/// HTML report generation via MiniJinja templates.
///
/// Templates are embedded at compile time (Jinja2-compatible syntax).
/// For the GUI: the rendered HTML is returned as a String and displayed
/// in a new Tauri WebView window. For the CLI: written to the results directory.
use crate::Result;
use matchy_core::SimulationResult;
use minijinja::{context, Environment};
use serde_json::Value;

// Embed templates at compile time
const REPORT_TEMPLATE: &str = include_str!("../templates/report.html");
const TRACE_REPORT_TEMPLATE: &str = include_str!("../templates/trace_report.html");

/// Render the standard simulation report to an HTML string.
///
/// `pedigree_image_b64`: base64-encoded PNG of the pedigree graph (from react-flow export)
/// `chart_images_b64`: map of chart_name → base64-encoded PNG (convergence plots)
pub fn render_report(
    result: &SimulationResult,
    pedigree_image_b64: Option<&str>,
    chart_images_b64: &std::collections::HashMap<String, String>,
) -> Result<String> {
    let mut env = Environment::new();
    env.add_template("report.html", REPORT_TEMPLATE)?;
    let tmpl = env.get_template("report.html")?;

    // Build template context from SimulationResult
    let ctx = context! {
        simulation_name => &result.parameters.simulation_name,
        user_name => &result.parameters.user_name,
        suspect => &result.parameters.suspect,
        trace_mode => result.parameters.trace_mode,
        converged => result.converged,
        trials => result.trials,
        pedigree_image => pedigree_image_b64.unwrap_or(""),
        inside_probabilities => result.inside_match_probabilities.as_ref().map(|p| {
            serde_json::to_value(p).unwrap_or(Value::Null)
        }),
        outside_probability => result.outside_match_probability.map(|p| p.to_string()),
        per_individual_probabilities => result.per_individual_probabilities.as_ref().map(|m| {
            serde_json::to_value(m).unwrap_or(Value::Null)
        }),
        charts => chart_images_b64,
    };

    Ok(tmpl.render(ctx)?)
}

/// Render a trace mode report to an HTML string.
pub fn render_trace_report(
    result: &SimulationResult,
    pedigree_image_b64: Option<&str>,
    chart_images_b64: &std::collections::HashMap<String, String>,
) -> Result<String> {
    let mut env = Environment::new();
    env.add_template("trace_report.html", TRACE_REPORT_TEMPLATE)?;
    let tmpl = env.get_template("trace_report.html")?;

    let ctx = context! {
        simulation_name => &result.parameters.simulation_name,
        user_name => &result.parameters.user_name,
        converged => result.converged,
        trials => result.trials,
        pedigree_image => pedigree_image_b64.unwrap_or(""),
        per_individual_probabilities => result.per_individual_probabilities.as_ref().map(|m| {
            // Sort by probability descending for display
            let mut entries: Vec<_> = m.iter().map(|(k, v)| (k.clone(), v.to_string())).collect();
            entries.sort_by(|a, b| b.1.partial_cmp(&a.1).unwrap_or(std::cmp::Ordering::Equal));
            entries
        }),
        charts => chart_images_b64,
    };

    Ok(tmpl.render(ctx)?)
}

/// Normalise per-individual probabilities to sum to 1 (trace mode).
pub fn normalize_probabilities(
    probs: &std::collections::HashMap<String, rust_decimal::Decimal>,
) -> std::collections::HashMap<String, rust_decimal::Decimal> {
    use rust_decimal::Decimal;
    let total: Decimal = probs.values().sum();
    if total == Decimal::ZERO {
        return probs.clone();
    }
    probs.iter().map(|(k, v)| (k.clone(), v / total)).collect()
}
