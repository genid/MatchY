import {
  forwardRef,
  useEffect,
  useImperativeHandle,
  useRef,
} from "react";
import { useT } from "../i18n";
import {
  Chart as ChartJS,
  CategoryScale,
  LinearScale,
  PointElement,
  LineElement,
  Title,
  Tooltip,
  Legend,
  type ChartData,
  type ChartOptions,
} from "chart.js";
import zoomPlugin from "chartjs-plugin-zoom";
import { Line } from "react-chartjs-2";
import type { ProgressEvent } from "../types/matchy";

ChartJS.register(
  CategoryScale,
  LinearScale,
  PointElement,
  LineElement,
  Title,
  Tooltip,
  Legend,
  zoomPlugin
);

interface Props {
  events: ProgressEvent[];
  stage: "pedigree_probability" | "extended_pedigree_probability" | "inside_match_probability" | "outside_match_probability";
  title: string;
  convergenceCriterion?: number;
  batchLength?: number;
}

export interface ConvergenceChartRef {
  toBase64Image: () => string | null;
}

const MODEL_COLORS = ["#3b82f6", "#10b981", "#f59e0b"];

export const ConvergenceChart = forwardRef<ConvergenceChartRef, Props>(
  ({ events, stage, title, convergenceCriterion, batchLength = 1 }, ref) => {
    const chartRef = useRef<ChartJS<"line"> | null>(null);
    const userZoomRange = useRef<{ min: number; max: number } | null>(null);
    const isProgrammatic = useRef(false);
    const t = useT();
    const MODEL_LABELS = [0, 1, 2].map((i) => `${t("conv_model_prefix")} ${i}`);

    useImperativeHandle(ref, () => ({
      toBase64Image: () => {
        if (!chartRef.current) return null;
        return chartRef.current.toBase64Image();
      },
    }));

    // Filter events for this stage, grouped by model
    const filteredByModel: ProgressEvent[][] = [[], [], []];
    for (const ev of events) {
      if (ev.stage === stage && ev.model >= 0 && ev.model <= 2) {
        filteredByModel[ev.model].push(ev);
      }
    }

    // X-axis: iteration counts (batchLength, 2×batchLength, ...) based on the longest model series
    const maxLen = Math.max(...filteredByModel.map((m) => m.length), 0);
    const labels = Array.from({ length: maxLen }, (_, i) => String((i + 1) * batchLength));

    // Latest values per model
    const latestValues = filteredByModel
      .map((m) => m.length > 0 ? parseFloat(m[m.length - 1].currentMean) : null)
      .filter((v): v is number => v !== null);

    const currentAvg = latestValues.length > 0
      ? latestValues.reduce((a, b) => a + b, 0) / latestValues.length
      : null;

    const currentVariance = latestValues.length > 1 && currentAvg && currentAvg > 0
      ? (Math.max(...latestValues) - Math.min(...latestValues)) / currentAvg
      : null;

    // Auto-zoom x-axis to show the latest 75% of data on each update.
    // Once the user pans or zooms manually their range is stored and restored
    // after every Chart.js data update. Reset Zoom clears the stored range.
    useEffect(() => {
      const chart = chartRef.current;
      if (!chart || maxLen === 0) return;

      isProgrammatic.current = true;
      if (userZoomRange.current !== null) {
        // Restore the user's manually chosen range after Chart.js re-renders.
        const { min, max } = userZoomRange.current;
        chart.zoomScale("x", { min, max }, "none");
      } else if (maxLen >= 4) {
        // Auto-zoom: skip the first 25% (likely still stochastic), show the rest.
        const startIdx = Math.floor(maxLen * 0.25);
        chart.zoomScale("x", { min: startIdx, max: maxLen - 1 }, "none");
      }
      isProgrammatic.current = false;
    });

    const modelDatasets = filteredByModel.map((modelEvents, modelIdx) => ({
      label: MODEL_LABELS[modelIdx],
      data: modelEvents.map((ev) => parseFloat(ev.currentMean)),
      borderColor: MODEL_COLORS[modelIdx],
      backgroundColor: MODEL_COLORS[modelIdx] + "33",
      borderWidth: 2,
      pointRadius: 2,
      tension: 0.2,
    }));

    // Reference line datasets (only when we have data)
    const refDatasets: ChartData<"line">["datasets"] = [];
    if (currentAvg !== null && maxLen > 0) {
      const avgArr = Array(maxLen).fill(currentAvg);
      refDatasets.push({
        label: t("conv_current_avg"),
        data: avgArr,
        borderColor: "#6366f1",
        backgroundColor: "transparent",
        borderWidth: 1.5,
        borderDash: [6, 4],
        pointRadius: 0,
        tension: 0,
      } as ChartData<"line">["datasets"][0]);

      if (convergenceCriterion !== undefined && convergenceCriterion > 0) {
        const upper = Array(maxLen).fill(currentAvg * (1 + convergenceCriterion));
        const lower = Array(maxLen).fill(currentAvg * (1 - convergenceCriterion));
        refDatasets.push({
          label: t("conv_criterion_upper").replace("{pct}", (convergenceCriterion * 100).toFixed(1)),
          data: upper,
          borderColor: "#ef4444",
          backgroundColor: "transparent",
          borderWidth: 1,
          borderDash: [3, 4],
          pointRadius: 0,
          tension: 0,
        } as ChartData<"line">["datasets"][0]);
        refDatasets.push({
          label: t("conv_criterion_lower"),
          data: lower,
          borderColor: "#ef4444",
          backgroundColor: "transparent",
          borderWidth: 1,
          borderDash: [3, 4],
          pointRadius: 0,
          tension: 0,
        } as ChartData<"line">["datasets"][0]);
      }
    }

    const data: ChartData<"line"> = { labels, datasets: [...modelDatasets, ...refDatasets] };

    const options: ChartOptions<"line"> = {
      responsive: true,
      maintainAspectRatio: false,
      animation: false,
      plugins: {
        legend: { position: "top" as const },
        title: {
          display: true,
          text: title,
          font: { size: 13 },
        },
        zoom: {
          pan: {
            enabled: true,
            mode: "xy",
            onPan: (ctx: { chart: ChartJS }) => {
              if (!isProgrammatic.current) {
                const xScale = ctx.chart.scales["x"];
                if (xScale) userZoomRange.current = { min: xScale.min, max: xScale.max };
              }
            },
          },
          zoom: {
            wheel: { enabled: true },
            pinch: { enabled: true },
            mode: "xy",
            onZoom: (ctx: { chart: ChartJS }) => {
              if (!isProgrammatic.current) {
                const xScale = ctx.chart.scales["x"];
                if (xScale) userZoomRange.current = { min: xScale.min, max: xScale.max };
              }
            },
          },
        },
      },
      scales: {
        x: {
          title: { display: true, text: t("conv_batch") },
          ticks: { maxTicksLimit: 10 },
        },
        y: {
          title: { display: true, text: t("conv_mean_probability") },
          ticks: {
            callback: (value) =>
              typeof value === "number" ? value.toExponential(2) : value,
          },
        },
      },
    };

    if (maxLen === 0) {
      return null;
    }

    // Use backend-reported converged flag as authoritative (Rust checks all running means,
    // frontend can only check final values which may differ).
    const convergedFromEvents = filteredByModel.some(
      (m) => m.length > 0 && m[m.length - 1].converged
    );
    const converged = convergedFromEvents
      ? true
      : convergenceCriterion !== undefined && currentVariance !== null
      ? currentVariance <= convergenceCriterion
      : null;

    return (
      <div>
        {/* Variance badge */}
        {currentVariance !== null && (
          <div className="flex items-center gap-3 mb-1 text-xs">
            <span className="text-gray-500">
              {t("conv_variance_label")}{" "}
              <span className={`font-mono font-semibold ${converged ? "text-emerald-600" : "text-orange-600"}`}>
                {(currentVariance * 100).toFixed(3)}%
              </span>
            </span>
            {convergenceCriterion !== undefined && (
              <span className="text-gray-400">
                {t("conv_criterion_label")} {(convergenceCriterion * 100).toFixed(1)}%
              </span>
            )}
            {converged !== null && (
              <span className={`font-medium ${converged ? "text-emerald-600" : "text-orange-500"}`}>
                {converged ? t("conv_converged") : t("conv_converging")}
              </span>
            )}
          </div>
        )}
        <div className="relative h-72">
          <button
            onClick={() => {
              userZoomRange.current = null;
              chartRef.current?.resetZoom();
            }}
            className="absolute top-1 right-1 z-10 text-xs bg-white border border-gray-300 hover:bg-gray-50 text-gray-600 px-2 py-0.5 rounded shadow-sm"
          >
            {t("conv_reset_zoom")}
          </button>
          <Line
            ref={chartRef}
            data={data}
            options={options}
            style={{ height: "100%" }}
          />
        </div>
      </div>
    );
  }
);

ConvergenceChart.displayName = "ConvergenceChart";
