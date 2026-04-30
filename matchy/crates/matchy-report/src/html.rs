use crate::Result;
use matchy_core::SimulationResult;
use minijinja::{context, Environment};
use rust_decimal::Decimal;
use serde_json::Value;
use std::collections::HashMap;

const REPORT_TEMPLATE: &str = include_str!("../templates/report.html");
const TRACE_REPORT_TEMPLATE: &str = include_str!("../templates/trace_report.html");
const LOGO_BYTES: &[u8] = include_bytes!("../assets/logo.png");

fn logo_data_url() -> String {
    format!("data:image/png;base64,{}", base64_encode(LOGO_BYTES))
}

fn base64_encode(data: &[u8]) -> String {
    const CHARS: &[u8] = b"ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz0123456789+/";
    let mut out = String::with_capacity((data.len() + 2) / 3 * 4);
    for chunk in data.chunks(3) {
        let b0 = chunk[0] as usize;
        let b1 = if chunk.len() > 1 { chunk[1] as usize } else { 0 };
        let b2 = if chunk.len() > 2 { chunk[2] as usize } else { 0 };
        out.push(CHARS[b0 >> 2] as char);
        out.push(CHARS[((b0 & 3) << 4) | (b1 >> 4)] as char);
        out.push(if chunk.len() > 1 { CHARS[((b1 & 0xf) << 2) | (b2 >> 6)] as char } else { '=' });
        out.push(if chunk.len() > 2 { CHARS[b2 & 0x3f] as char } else { '=' });
    }
    out
}

// ---------------------------------------------------------------------------
// Helpers
// ---------------------------------------------------------------------------

/// 3 significant figures: fixed notation for [1e-3, 1e3), scientific otherwise.
fn fmt_3sig(f: f64) -> String {
    if f == 0.0 { return "0".to_string(); }
    let abs = f.abs();
    if abs >= 1e-3 && abs < 1e3 {
        let mag = abs.log10().floor() as i32;
        let decimals = (2 - mag).max(0) as usize;
        format!("{:.prec$}", f, prec = decimals)
    } else {
        format!("{:.2E}", f)
    }
}

fn fmt_decimal(d: Decimal) -> String {
    fmt_3sig(f64::try_from(d).unwrap_or(0.0))
}

fn fmt_pct(d: Decimal) -> String {
    fmt_3sig(f64::try_from(d).unwrap_or(0.0) * 100.0)
}

fn fmt_lr(prob: Decimal) -> String {
    let f = f64::try_from(prob).unwrap_or(0.0);
    if f <= 0.0 { return "∞".to_string(); }
    fmt_3sig(1.0 / f)
}

/// Sort a HashMap<String, Decimal> descending by value → [(name, prob_str, pct_str, lr_str)]
fn sorted_per_individual(
    map: &HashMap<String, Decimal>,
) -> (Vec<(String, String, String, String)>, String) {
    let mut entries: Vec<(String, Decimal)> =
        map.iter().map(|(k, v)| (k.clone(), *v)).collect();
    entries.sort_by(|a, b| b.1.cmp(&a.1));

    let mut prob_sum = 0.0f64;
    let mut prob_count = 0usize;

    let rows = entries
        .into_iter()
        .map(|(name, prob)| {
            let lr_str = fmt_lr(prob);
            let f = f64::try_from(prob).unwrap_or(0.0);
            if f > 0.0 {
                prob_sum += f;
                prob_count += 1;
            }
            (name, fmt_decimal(prob), fmt_pct(prob), lr_str)
        })
        .collect();

    let avg_lr = if prob_count > 0 && prob_sum > 0.0 {
        fmt_3sig(prob_count as f64 / prob_sum)
    } else {
        String::new()
    };

    (rows, avg_lr)
}

/// Parse pedigree_json → list of {name, haplotype_class, exclude}
fn parse_individuals(pedigree_json: Option<&str>) -> Vec<Value> {
    let json = match pedigree_json.and_then(|s| serde_json::from_str::<Value>(s).ok()) {
        Some(v) => v,
        None => return vec![],
    };
    json["individuals"]
        .as_array()
        .cloned()
        .unwrap_or_default()
        .into_iter()
        .filter_map(|ind| {
            let name = ind["name"].as_str()?.to_string();
            let hc = ind["haplotypeClass"].as_str().unwrap_or("unknown").to_string();
            let exclude = ind["exclude"].as_bool().unwrap_or(false);
            Some(serde_json::json!({
                "name": name,
                "haplotype_class": hc,
                "exclude": if exclude { "yes" } else { "no" },
            }))
        })
        .collect()
}

/// Return individual names from haplotypeTable keys (sorted, TRACE excluded).
/// This is the correct source because the pedigree JSON in the store has stale
/// haplotypeClass="unknown" for everyone after haplotype loading.
fn haplotype_table_names(haplotypes_json: Option<&str>) -> Vec<String> {
    let val = match haplotypes_json.and_then(|s| serde_json::from_str::<Value>(s).ok()) {
        Some(v) => v,
        None => return vec![],
    };
    let mut names: Vec<String> = val["haplotypeTable"]
        .as_object()
        .map(|table| {
            table
                .keys()
                .filter(|k| *k != "TRACE")
                .cloned()
                .collect()
        })
        .unwrap_or_default();
    names.sort();
    names
}


/// Build haplotype table rows: [{marker, alleles:[str], highlight:bool}]
/// alleles are in same order as `known_names`.
fn build_haplotype_rows(
    haplotypes_json: Option<&str>,
    markers_json: Option<&str>,
    known_names: &[String],
    trace_col: bool,
) -> Vec<Value> {
    let hap_val = match haplotypes_json.and_then(|s| serde_json::from_str::<Value>(s).ok()) {
        Some(v) => v,
        None => return vec![],
    };
    let marker_names = if let Some(ms) = markers_json {
        serde_json::from_str::<Value>(ms)
            .ok()
            .and_then(|v| {
                v.as_array().map(|arr| {
                    arr.iter()
                        .filter_map(|m| m["name"].as_str().map(|s| s.to_string()))
                        .collect::<Vec<_>>()
                })
            })
            .unwrap_or_default()
    } else {
        // Fallback: use markerNames from haplotypes JSON
        hap_val["markerNames"]
            .as_array()
            .cloned()
            .unwrap_or_default()
            .into_iter()
            .filter_map(|v| v.as_str().map(|s| s.to_string()))
            .collect()
    };

    let table = &hap_val["haplotypeTable"];
    let trace = &hap_val["traceHaplotype"];

    marker_names
        .iter()
        .map(|marker| {
            let trace_allele = if trace_col {
                trace[marker].as_str().unwrap_or("—").to_string()
            } else {
                String::new()
            };

            let alleles: Vec<String> = known_names
                .iter()
                .map(|name| {
                    table[name][marker]
                        .as_str()
                        .map(|s| s.to_string())
                        .unwrap_or_else(|| "—".to_string())
                })
                .collect();

            // Highlight if any non-dash allele differs
            let non_dash: Vec<&str> = if trace_col {
                let mut v: Vec<&str> = alleles
                    .iter()
                    .filter(|a| *a != "—")
                    .map(|a| a.as_str())
                    .collect();
                if trace_allele != "—" && !trace_allele.is_empty() {
                    v.push(&trace_allele);
                }
                v
            } else {
                alleles
                    .iter()
                    .filter(|a| *a != "—")
                    .map(|a| a.as_str())
                    .collect()
            };
            let unique: std::collections::HashSet<&str> = non_dash.iter().copied().collect();
            let highlight = unique.len() > 1;

            serde_json::json!({
                "marker": marker,
                "trace_allele": trace_allele,
                "alleles": alleles,
                "highlight": highlight,
            })
        })
        .collect()
}

/// Build marker rows: [{name, rate, copies}]
fn build_marker_rows(markers_json: Option<&str>) -> Vec<Value> {
    let json = match markers_json.and_then(|s| serde_json::from_str::<Value>(s).ok()) {
        Some(v) => v,
        None => return vec![],
    };
    let mut rows: Vec<Value> = json
        .as_array()
        .cloned()
        .unwrap_or_default()
        .into_iter()
        .filter_map(|m| {
            let name = m["name"].as_str()?.to_string();
            let rate = m["mutationRate"].as_f64().unwrap_or(0.0);
            let copies = m["numberOfCopies"].as_u64().unwrap_or(1);
            Some(serde_json::json!({
                "name": name,
                "rate": fmt_3sig(rate),
                "copies": copies,
            }))
        })
        .collect();
    rows.sort_by(|a, b| {
        a["name"]
            .as_str()
            .unwrap_or("")
            .cmp(b["name"].as_str().unwrap_or(""))
    });
    rows
}

/// Build chart data JSON string for interactive Chart.js charts in the report.
/// Input is the serialised ProgressEvent[] from the frontend (camelCase keys).
fn build_chart_data_json(events_json: Option<&str>) -> String {
    let events: Vec<Value> = match events_json
        .and_then(|s| serde_json::from_str::<Value>(s).ok())
    {
        Some(Value::Array(v)) => v,
        _ => return "null".to_string(),
    };

    let stages = [
        "pedigree_probability",
        "extended_pedigree_probability",
        "inside_match_probability",
        "outside_match_probability",
    ];

    let mut result = serde_json::Map::new();
    for stage in &stages {
        let mut models: [Vec<f64>; 3] = [vec![], vec![], vec![]];
        for ev in &events {
            if ev["stage"].as_str().unwrap_or("") != *stage {
                continue;
            }
            let model = ev["model"].as_u64().unwrap_or(0) as usize;
            let mean = ev["currentMean"]
                .as_str()
                .unwrap_or("0")
                .parse::<f64>()
                .unwrap_or(0.0);
            if model < 3 {
                models[model].push(mean);
            }
        }
        if models.iter().any(|m| !m.is_empty()) {
            result.insert(
                stage.to_string(),
                serde_json::json!({
                    "0": models[0],
                    "1": models[1],
                    "2": models[2],
                }),
            );
        }
    }

    if result.is_empty() {
        "null".to_string()
    } else {
        serde_json::to_string(&Value::Object(result)).unwrap_or_else(|_| "null".to_string())
    }
}

fn sorted_charts(chart_images_b64: &HashMap<String, String>) -> Value {
    let mut list: Vec<(String, String)> = chart_images_b64
        .iter()
        .map(|(k, v)| (k.clone(), v.clone()))
        .collect();
    list.sort_by(|a, b| a.0.cmp(&b.0));
    serde_json::to_value(list).unwrap_or(Value::Null)
}

// ---------------------------------------------------------------------------
// render_report
// ---------------------------------------------------------------------------

pub fn render_report(
    result: &SimulationResult,
    pedigree_image_b64: Option<&str>,
    chart_images_b64: &HashMap<String, String>,
    pedigree_json: Option<&str>,
    haplotypes_json: Option<&str>,
    markers_json: Option<&str>,
    report_date: Option<&str>,
    progress_events_json: Option<&str>,
    lang_json: Option<&str>,
) -> Result<String> {
    let mut env = Environment::new();
    env.add_template("report.html", REPORT_TEMPLATE)?;
    let tmpl = env.get_template("report.html")?;

    // Inside match probability — P(at least one other pedigree member matches) = sum of all k
    let (inside_prob, inside_k1_pct, inside_k1_lr) =
        if let Some(ref p) = result.inside_match_probabilities {
            let total: Decimal = p.probabilities.values().sum();
            if total > Decimal::ZERO {
                let f = f64::try_from(total).unwrap_or(0.0);
                let pct = fmt_3sig(f * 100.0);
                let lr = if f > 0.0 { fmt_3sig((1.0 - f) / f) } else { "∞".to_string() };
                (Some(fmt_decimal(total)), Some(pct), Some(lr))
            } else {
                (None, None, None)
            }
        } else {
            (None, None, None)
        };

    // --- Outside probability ---
    let outside_prob: Value = result
        .outside_match_probability
        .map(|p| Value::String(fmt_decimal(p)))
        .unwrap_or(Value::Null);

    let outside_lr: Value = result
        .outside_match_probability
        .map(|p| {
            let f = f64::try_from(p).unwrap_or(0.0);
            if f > 0.0 { Value::String(fmt_3sig(1.0 / f)) } else { Value::String("∞".to_string()) }
        })
        .unwrap_or(Value::Null);

    // --- Per-individual ---
    let (per_individual_rows, per_individual_avg_lr) =
        result.per_individual_probabilities.as_ref()
            .map(|m| sorted_per_individual(m))
            .unwrap_or_default();
    let per_individual: Value = if per_individual_rows.is_empty() {
        Value::Null
    } else {
        serde_json::to_value(&per_individual_rows).unwrap_or(Value::Null)
    };
    let per_individual_avg_lr = if per_individual_avg_lr.is_empty() {
        Value::Null
    } else {
        Value::String(per_individual_avg_lr)
    };

    // --- Pedigree info ---
    let pedigree_individuals =
        serde_json::to_value(parse_individuals(pedigree_json)).unwrap_or(Value::Null);

    // --- Haplotype table ---
    // Derive known individual names from the haplotypeTable keys (individuals that actually
    // have haplotypes), sorted. The pedigree JSON has stale haplotypeClass="unknown" for
    // everyone because the store isn't updated after haplotype loading.
    let known_names = haplotype_table_names(haplotypes_json);
    let haplotype_rows = serde_json::to_value(
        build_haplotype_rows(haplotypes_json, markers_json, &known_names, false)
    ).unwrap_or(Value::Null);
    let known_ind_names = serde_json::to_value(&known_names).unwrap_or(Value::Null);

    // --- Marker table ---
    let marker_rows = serde_json::to_value(build_marker_rows(markers_json)).unwrap_or(Value::Null);

    // --- Params ---
    let bias_str = if result.parameters.adaptive_bias {
        "Adaptive".to_string()
    } else {
        result.parameters.bias
            .map(|b| format!("{:.2}", b))
            .unwrap_or_else(|| "Auto".to_string())
    };

    // Separate extended pedigree from the convergence chart images.
    let extended_pedigree_image = chart_images_b64.get("extended_pedigree").cloned().unwrap_or_default();
    let convergence_charts: HashMap<String, String> = chart_images_b64
        .iter()
        .filter(|(k, _)| *k != "extended_pedigree")
        .map(|(k, v)| (k.clone(), v.clone()))
        .collect();

    // --- Performance stats ---
    let pedigree_iters_per_model = result.pedigree_stats.iterations_per_model.clone();
    let pedigree_runtime = result.pedigree_stats.formatted_runtime();
    let pedigree_model_probs = result.pedigree_stats.model_probabilities.clone();

    let inside_iters_per_model = result.inside_stats.as_ref().map(|s| s.iterations_per_model.clone()).unwrap_or_default();
    let inside_runtime = result.inside_stats.as_ref().map(|s| s.formatted_runtime()).unwrap_or_default();
    let inside_model_probs = result.inside_stats.as_ref().map(|s| s.model_probabilities.clone()).unwrap_or_default();

    let ext_ped_iters_per_model = result.extended_pedigree_stats.as_ref().map(|s| s.iterations_per_model.clone()).unwrap_or_default();
    let ext_ped_runtime = result.extended_pedigree_stats.as_ref().map(|s| s.formatted_runtime()).unwrap_or_default();
    let ext_ped_model_probs = result.extended_pedigree_stats.as_ref().map(|s| s.model_probabilities.clone()).unwrap_or_default();

    let outside_iters_per_model = result.outside_stats.as_ref().map(|s| s.iterations_per_model.clone()).unwrap_or_default();
    let outside_runtime = result.outside_stats.as_ref().map(|s| s.formatted_runtime()).unwrap_or_default();
    let outside_model_probs = result.outside_stats.as_ref().map(|s| s.model_probabilities.clone()).unwrap_or_default();

    let total_iterations: u64 = result.pedigree_stats.total_iterations()
        + result.inside_stats.as_ref().map(|s| s.total_iterations()).unwrap_or(0)
        + result.extended_pedigree_stats.as_ref().map(|s| s.total_iterations()).unwrap_or(0)
        + result.outside_stats.as_ref().map(|s| s.total_iterations()).unwrap_or(0);

    let total_runtime = {
        let s = result.total_runtime_secs;
        if s < 60.0 { format!("{:.1}s", s) } else { format!("{:.0}m {:.0}s", (s / 60.0).floor(), s % 60.0) }
    };

    let avg_pedigree_prob = result.inside_match_probabilities
        .as_ref()
        .map(|p| fmt_decimal(p.average_pedigree_probability))
        .unwrap_or_default();

    // Consecutive footnote numbers — only assigned when the corresponding row is shown.
    let mut fn_n: usize = 0;
    let fn1: usize = if !avg_pedigree_prob.is_empty() { fn_n += 1; fn_n } else { 0 };
    let fn2: usize = if inside_prob.is_some() { fn_n += 1; fn_n } else { 0 };
    let fn3: usize = if !ext_ped_iters_per_model.is_empty() { fn_n += 1; fn_n } else { 0 };
    let fn4: usize = if result.outside_match_probability.is_some() { fn_n += 1; fn_n } else { 0 };

    let ctx = context! {
        simulation_name => &result.parameters.simulation_name,
        user_name => &result.parameters.user_name,
        suspect => &result.parameters.suspect,
        converged => result.converged,
        trials => result.trials,
        date => report_date.unwrap_or(""),
        avg_pedigree_prob => avg_pedigree_prob,
        inside_prob => inside_prob,
        inside_k1_pct => inside_k1_pct,
        inside_k1_lr => inside_k1_lr,
        outside_prob => outside_prob,
        outside_lr => outside_lr,
        per_individual => per_individual,
        per_individual_avg_lr => per_individual_avg_lr,
        pedigree_image => pedigree_image_b64.unwrap_or(""),
        extended_pedigree_image => extended_pedigree_image,
        pedigree_individuals => pedigree_individuals,
        known_ind_names => known_ind_names,
        haplotype_rows => haplotype_rows,
        marker_rows => marker_rows,
        params_two_step => format!("{:.4}", result.parameters.two_step_mutation_fraction),
        params_batch => result.parameters.batch_length,
        params_convergence => format!("{:.4}", result.parameters.convergence_criterion),
        params_bias => bias_str,
        charts => sorted_charts(&convergence_charts),
        logo => logo_data_url(),
        chart_data_json => build_chart_data_json(progress_events_json),
        pedigree_iters_per_model => pedigree_iters_per_model,
        pedigree_runtime => pedigree_runtime,
        pedigree_model_probs => pedigree_model_probs,
        inside_iters_per_model => inside_iters_per_model,
        inside_runtime => inside_runtime,
        inside_model_probs => inside_model_probs,
        ext_ped_iters_per_model => ext_ped_iters_per_model,
        ext_ped_runtime => ext_ped_runtime,
        ext_ped_model_probs => ext_ped_model_probs,
        outside_iters_per_model => outside_iters_per_model,
        outside_runtime => outside_runtime,
        outside_model_probs => outside_model_probs,
        total_iterations => total_iterations,
        total_runtime => total_runtime,
        fn1 => fn1,
        fn2 => fn2,
        fn3 => fn3,
        fn4 => fn4,
        lang => lang_json.and_then(|s| serde_json::from_str::<Value>(s).ok()).unwrap_or(Value::Null),
    };

    Ok(tmpl.render(ctx)?)
}

// ---------------------------------------------------------------------------
// render_trace_report
// ---------------------------------------------------------------------------

pub fn render_trace_report(
    result: &SimulationResult,
    pedigree_image_b64: Option<&str>,
    chart_images_b64: &HashMap<String, String>,
    pedigree_json: Option<&str>,
    haplotypes_json: Option<&str>,
    markers_json: Option<&str>,
    report_date: Option<&str>,
    progress_events_json: Option<&str>,
    lang_json: Option<&str>,
) -> Result<String> {
    let mut env = Environment::new();
    env.add_template("trace_report.html", TRACE_REPORT_TEMPLATE)?;
    let tmpl = env.get_template("trace_report.html")?;

    // --- Build ranked list ---
    let mut ranked: Vec<(String, Decimal)> = Vec::new();
    if let Some(ref m) = result.per_individual_probabilities {
        for (name, prob) in m {
            ranked.push((name.clone(), *prob));
        }
    }
    if let Some(outside) = result.outside_match_probability {
        ranked.push(("Outside Pedigree".to_string(), outside));
    }
    ranked.sort_by(|a, b| b.1.cmp(&a.1));

    let max_prob = ranked.first().map(|(_, v)| *v).unwrap_or(Decimal::ZERO);
    let most_likely_donor = ranked.first().map(|(n, _)| n.clone()).unwrap_or_default();

    let ranked_individuals: Vec<(String, String, String)> = ranked
        .iter()
        .map(|(name, prob)| {
            let pct = if max_prob > Decimal::ZERO {
                let ratio = prob / max_prob;
                format!("{:.1}", f64::try_from(ratio * Decimal::from(100)).unwrap_or(0.0))
            } else {
                "0.0".to_string()
            };
            (name.clone(), fmt_decimal(*prob), pct)
        })
        .collect();

    // --- Pedigree info ---
    let pedigree_individuals =
        serde_json::to_value(parse_individuals(pedigree_json)).unwrap_or(Value::Null);

    // --- Haplotype table (with TRACE column) ---
    let has_trace = result.parameters.trace_mode;
    let known_names = haplotype_table_names(haplotypes_json);
    let haplotype_rows = serde_json::to_value(
        build_haplotype_rows(haplotypes_json, markers_json, &known_names, has_trace)
    ).unwrap_or(Value::Null);
    let known_ind_names = serde_json::to_value(&known_names).unwrap_or(Value::Null);

    // --- Marker table ---
    let marker_rows = serde_json::to_value(build_marker_rows(markers_json)).unwrap_or(Value::Null);

    // --- Params ---
    let bias_str = if result.parameters.adaptive_bias {
        "Adaptive".to_string()
    } else {
        result.parameters.bias
            .map(|b| format!("{:.2}", b))
            .unwrap_or_else(|| "Auto".to_string())
    };

    let extended_pedigree_image = chart_images_b64.get("extended_pedigree").cloned().unwrap_or_default();
    let convergence_charts: HashMap<String, String> = chart_images_b64
        .iter()
        .filter(|(k, _)| *k != "extended_pedigree")
        .map(|(k, v)| (k.clone(), v.clone()))
        .collect();

    // --- Performance stats ---
    let trace_pedigree_iters_per_model = result.pedigree_stats.iterations_per_model.clone();
    let trace_pedigree_runtime = result.pedigree_stats.formatted_runtime();
    let trace_pedigree_model_probs = result.pedigree_stats.model_probabilities.clone();

    let trace_inside_iters_per_model = result.inside_stats.as_ref().map(|s| s.iterations_per_model.clone()).unwrap_or_default();
    let trace_inside_runtime = result.inside_stats.as_ref().map(|s| s.formatted_runtime()).unwrap_or_default();
    let trace_inside_model_probs = result.inside_stats.as_ref().map(|s| s.model_probabilities.clone()).unwrap_or_default();

    let trace_ext_ped_iters_per_model = result.extended_pedigree_stats.as_ref().map(|s| s.iterations_per_model.clone()).unwrap_or_default();
    let trace_ext_ped_runtime = result.extended_pedigree_stats.as_ref().map(|s| s.formatted_runtime()).unwrap_or_default();
    let trace_ext_ped_model_probs = result.extended_pedigree_stats.as_ref().map(|s| s.model_probabilities.clone()).unwrap_or_default();

    let trace_outside_iters_per_model = result.outside_stats.as_ref().map(|s| s.iterations_per_model.clone()).unwrap_or_default();
    let trace_outside_runtime = result.outside_stats.as_ref().map(|s| s.formatted_runtime()).unwrap_or_default();
    let trace_outside_model_probs = result.outside_stats.as_ref().map(|s| s.model_probabilities.clone()).unwrap_or_default();

    let trace_total_iterations: u64 = result.pedigree_stats.total_iterations()
        + result.inside_stats.as_ref().map(|s| s.total_iterations()).unwrap_or(0)
        + result.extended_pedigree_stats.as_ref().map(|s| s.total_iterations()).unwrap_or(0)
        + result.outside_stats.as_ref().map(|s| s.total_iterations()).unwrap_or(0);

    let trace_total_runtime = {
        let s = result.total_runtime_secs;
        if s < 60.0 { format!("{:.1}s", s) } else { format!("{:.0}m {:.0}s", (s / 60.0).floor(), s % 60.0) }
    };

    let ctx = context! {
        simulation_name => &result.parameters.simulation_name,
        user_name => &result.parameters.user_name,
        converged => result.converged,
        trials => result.trials,
        date => report_date.unwrap_or(""),
        most_likely_donor => most_likely_donor,
        ranked_individuals => serde_json::to_value(ranked_individuals).unwrap_or(Value::Null),
        has_trace => has_trace,
        pedigree_individuals => pedigree_individuals,
        known_ind_names => known_ind_names,
        haplotype_rows => haplotype_rows,
        marker_rows => marker_rows,
        params_two_step => format!("{:.4}", result.parameters.two_step_mutation_fraction),
        params_batch => result.parameters.batch_length,
        params_convergence => format!("{:.4}", result.parameters.convergence_criterion),
        params_bias => bias_str,
        pedigree_image => pedigree_image_b64.unwrap_or(""),
        extended_pedigree_image => extended_pedigree_image,
        charts => sorted_charts(&convergence_charts),
        logo => logo_data_url(),
        chart_data_json => build_chart_data_json(progress_events_json),
        pedigree_iters_per_model => trace_pedigree_iters_per_model,
        pedigree_runtime => trace_pedigree_runtime,
        pedigree_model_probs => trace_pedigree_model_probs,
        inside_iters_per_model => trace_inside_iters_per_model,
        inside_runtime => trace_inside_runtime,
        inside_model_probs => trace_inside_model_probs,
        ext_ped_iters_per_model => trace_ext_ped_iters_per_model,
        ext_ped_runtime => trace_ext_ped_runtime,
        ext_ped_model_probs => trace_ext_ped_model_probs,
        outside_iters_per_model => trace_outside_iters_per_model,
        outside_runtime => trace_outside_runtime,
        outside_model_probs => trace_outside_model_probs,
        total_iterations => trace_total_iterations,
        total_runtime => trace_total_runtime,
        lang => lang_json.and_then(|s| serde_json::from_str::<Value>(s).ok()).unwrap_or(Value::Null),
    };

    Ok(tmpl.render(ctx)?)
}

// ---------------------------------------------------------------------------
// normalize_probabilities (kept for potential external use)
// ---------------------------------------------------------------------------

pub fn normalize_probabilities(
    probs: &HashMap<String, rust_decimal::Decimal>,
) -> HashMap<String, rust_decimal::Decimal> {
    let total: Decimal = probs.values().sum();
    if total == Decimal::ZERO {
        return probs.clone();
    }
    probs.iter().map(|(k, v)| (k.clone(), v / total)).collect()
}
