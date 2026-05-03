## What's new in v2.0

v2.0 is a complete rewrite of MatchY in Rust, replacing the Python/Streamlit implementation (v1.0).

**Key improvements:**
- **Performance**: Significantly faster simulation engine with multi-threaded Monte Carlo sampling and importance weighting, reducing execution time by an order of magnitude on typical cases
- **Robustness**: Greater numerical stability and reliability across edge cases; more rigorous convergence detection
- **Graphical user interface**: New desktop application (Windows, Linux, macOS) built with Tauri — more intuitive and user-friendly than the previous browser-based Streamlit interface
- **Multi-language support**: The GUI is available in English, Dutch, German, French, Spanish, Portuguese, and Chinese
- **Cross-platform CLI**: A zero-dependency command-line binary for batch processing and scripted workflows

The underlying mathematical framework and algorithms are unchanged from v1.0: Monte Carlo simulation with importance sampling, the same mutation model (step changes with bias), and the same likelihood ratio framework for within- and outside-pedigree hypotheses.

The Python/Streamlit implementation is preserved in `python/` for reference.
