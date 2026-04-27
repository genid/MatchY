import { useEffect, useState } from "react";
import { useT } from "../i18n";

const TOUR_KEY = "matchy_tour_done";

const TOUR_ICONS = ["🧬", "🌳", "📊", "🔬", "⚙️", "📈", "💾"];

export function GuidedTour({ onClose }: { onClose: () => void }) {
  const [step, setStep] = useState(0);
  const t = useT();

  const steps = [
    { title: t("tour_step_0_title"), content: t("tour_step_0_content"), icon: TOUR_ICONS[0] },
    { title: t("tour_step_1_title"), content: t("tour_step_1_content"), icon: TOUR_ICONS[1] },
    { title: t("tour_step_2_title"), content: t("tour_step_2_content"), icon: TOUR_ICONS[2] },
    { title: t("tour_step_3_title"), content: t("tour_step_3_content"), icon: TOUR_ICONS[3] },
    { title: t("tour_step_4_title"), content: t("tour_step_4_content"), icon: TOUR_ICONS[4] },
    { title: t("tour_step_5_title"), content: t("tour_step_5_content"), icon: TOUR_ICONS[5] },
    { title: t("tour_step_6_title"), content: t("tour_step_6_content"), icon: TOUR_ICONS[6] },
  ];

  const current = steps[step];
  const isLast = step === steps.length - 1;

  const finish = () => {
    localStorage.setItem(TOUR_KEY, "1");
    onClose();
  };

  return (
    <div className="fixed inset-0 z-50 flex items-center justify-center bg-black/40">
      <div className="bg-white rounded-xl shadow-2xl w-full max-w-md mx-4 overflow-hidden">
        {/* Progress bar */}
        <div className="h-1 bg-gray-100">
          <div
            className="h-1 bg-blue-500 transition-all duration-300"
            style={{ width: `${((step + 1) / steps.length) * 100}%` }}
          />
        </div>

        <div className="p-6">
          <div className="text-4xl mb-3">{current.icon}</div>
          <h2 className="text-lg font-bold text-gray-900 mb-2">{current.title}</h2>
          <p className="text-sm text-gray-600 leading-relaxed">{current.content}</p>
        </div>

        <div className="px-6 pb-5 flex items-center justify-between">
          <button
            onClick={finish}
            className="text-xs text-gray-400 hover:text-gray-600 transition-colors"
          >
            {t("tour_skip")}
          </button>
          <div className="flex items-center gap-2">
            {/* Step dots */}
            <div className="flex gap-1 mr-3">
              {steps.map((_, i) => (
                <button
                  key={i}
                  onClick={() => setStep(i)}
                  className={`w-1.5 h-1.5 rounded-full transition-colors ${
                    i === step ? "bg-blue-500" : "bg-gray-300"
                  }`}
                />
              ))}
            </div>
            {step > 0 && (
              <button
                onClick={() => setStep((s) => s - 1)}
                className="text-sm border border-gray-300 text-gray-600 hover:bg-gray-50 px-3 py-1.5 rounded transition-colors"
              >
                {t("tour_back")}
              </button>
            )}
            <button
              onClick={() => (isLast ? finish() : setStep((s) => s + 1))}
              className="text-sm bg-blue-600 hover:bg-blue-700 text-white font-medium px-4 py-1.5 rounded transition-colors"
            >
              {isLast ? t("tour_get_started") : t("tour_next")}
            </button>
          </div>
        </div>
      </div>
    </div>
  );
}

export function useTour() {
  const [show, setShow] = useState(false);

  useEffect(() => {
    if (!localStorage.getItem(TOUR_KEY)) {
      setShow(true);
    }
  }, []);

  return {
    show,
    open: () => setShow(true),
    close: () => setShow(false),
  };
}
