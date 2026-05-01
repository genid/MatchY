use serde::Serialize;
use tauri::command;

#[derive(Serialize)]
#[serde(rename_all = "camelCase")]
pub struct UpdateInfo {
    pub latest_version: String,
    pub is_newer: bool,
    pub release_notes: String,
    pub download_url: String,
}

#[command]
pub async fn check_for_updates(current_version: String) -> Result<UpdateInfo, String> {
    let client = reqwest::Client::builder()
        .user_agent("MatchY-desktop")
        .build()
        .map_err(|e| e.to_string())?;

    let resp: serde_json::Value = client
        .get("https://api.github.com/repos/genid/MatchY/releases/latest")
        .send()
        .await
        .map_err(|e| e.to_string())?
        .json()
        .await
        .map_err(|e| e.to_string())?;

    let latest = resp["tag_name"]
        .as_str()
        .unwrap_or("")
        .trim_start_matches('v')
        .to_string();

    let notes = resp["body"].as_str().unwrap_or("").to_string();
    let url = resp["html_url"]
        .as_str()
        .unwrap_or("https://github.com/genid/MatchY/releases/latest")
        .to_string();

    let is_newer = is_version_newer(&latest, &current_version);

    Ok(UpdateInfo {
        latest_version: latest,
        is_newer,
        release_notes: notes,
        download_url: url,
    })
}

fn is_version_newer(latest: &str, current: &str) -> bool {
    let parse = |v: &str| -> (u32, u32, u32) {
        let parts: Vec<u32> = v.split('.').filter_map(|p| p.parse().ok()).collect();
        (
            parts.first().copied().unwrap_or(0),
            parts.get(1).copied().unwrap_or(0),
            parts.get(2).copied().unwrap_or(0),
        )
    };
    parse(latest) > parse(current)
}
