# Contributing to MatchY

Thank you for your interest in contributing to MatchY! We welcome contributions from the community.

## Table of Contents

- [Code of Conduct](#code-of-conduct)
- [Getting Started](#getting-started)
- [Development Setup](#development-setup)
- [How to Contribute](#how-to-contribute)
- [Coding Standards](#coding-standards)
- [Testing](#testing)
- [Pull Request Process](#pull-request-process)
- [Reporting Bugs](#reporting-bugs)
- [Suggesting Features](#suggesting-features)

## Code of Conduct

This project adheres to a Code of Conduct. By participating, you are expected to uphold this code. Please read [CODE_OF_CONDUCT.md](CODE_OF_CONDUCT.md) before contributing.

## Getting Started

1. Fork the repository on GitHub
2. Clone your fork locally:
   ```bash
   git clone https://github.com/your-username/MatchY.git
   cd MatchY
   ```
3. Add the upstream repository:
   ```bash
   git remote add upstream https://github.com/genid/MatchY.git
   ```

## Development Setup

### Prerequisites

- Python 3.9 or higher
- pip package manager
- Git

### Installation

1. Create a virtual environment:
   ```bash
   python -m venv venv
   source venv/bin/activate  # On Windows: venv\Scripts\activate
   ```

2. Install dependencies:
   ```bash
   pip install -r requirements.txt
   ```

3. Install development dependencies (optional):
   ```bash
   pip install black flake8 isort pytest
   ```

4. Install MatchY in editable mode:
   ```bash
   pip install -e .
   ```

### Running MatchY

- **GUI Mode:**
  ```bash
  streamlit run 1_🧬_Home.py
  ```

- **CLI Mode:**
  ```bash
  python main.py --config examples/config.ini
  ```

## How to Contribute

### Finding Issues to Work On

- Check the [Issues](https://github.com/genid/MatchY/issues) page
- Look for issues labeled `good first issue` or `help wanted`
- Feel free to ask questions in issue comments

### Creating a Branch

Create a descriptive branch name:

```bash
git checkout -b feature/add-new-marker-set
git checkout -b fix/haplotype-editor-bug
git checkout -b docs/update-user-manual
```

Branch naming conventions:
- `feature/` - New features
- `fix/` - Bug fixes
- `docs/` - Documentation changes
- `refactor/` - Code refactoring
- `test/` - Adding or updating tests

## Coding Standards

### Python Style Guide

We follow [PEP 8](https://www.python.org/dev/peps/pep-0008/) Python style guidelines.

**Key points:**
- Use 4 spaces for indentation (not tabs)
- Maximum line length: 100 characters
- Use descriptive variable and function names
- Add docstrings to functions and classes
- Use type hints where appropriate

### Code Formatting

We recommend using automated formatters:

```bash
# Format code with Black
black .

# Sort imports with isort
isort .

# Check style with flake8
flake8 .
```

### Documentation

- Update docstrings when modifying functions
- Update USER_MANUAL.md for user-facing changes
- Update PARAMETERS_REFERENCE.md when adding/modifying parameters
- Add inline comments for complex logic

## Testing

### Running Tests

```bash
# Run all tests
pytest

# Run specific test file
pytest tests/test_pedigree.py

# Run with coverage
pytest --cov=pedigree_lr
```

### Writing Tests

- Add tests for new features
- Ensure tests pass before submitting PR
- Aim for good code coverage
- Test edge cases and error conditions

## Pull Request Process

1. **Update your fork:**
   ```bash
   git fetch upstream
   git rebase upstream/master
   ```

2. **Make your changes:**
   - Write clear, concise commit messages
   - Follow the coding standards
   - Add tests if applicable
   - Update documentation

3. **Commit your changes:**
   ```bash
   git add .
   git commit -m "Add feature: description of changes"
   ```

4. **Push to your fork:**
   ```bash
   git push origin feature/your-feature-name
   ```

5. **Create a Pull Request:**
   - Go to the [MatchY repository](https://github.com/genid/MatchY)
   - Click "New Pull Request"
   - Select your fork and branch
   - Fill out the PR template with detailed description
   - Link related issues using keywords (e.g., "Fixes #123")

6. **PR Review Process:**
   - Maintainers will review your PR
   - Address any requested changes
   - Once approved, your PR will be merged

### Commit Message Guidelines

Write clear commit messages:

```
Add trace mode support for multiple markers

- Implement marker selection in trace analysis
- Update report template to show marker-specific results
- Add validation for marker compatibility

Fixes #123
```

Format:
- First line: Brief summary (50 chars or less)
- Blank line
- Detailed description (wrap at 72 chars)
- Reference issues at the end

## Reporting Bugs

Found a bug? Please create an issue using the bug report template:

1. Go to [Issues](https://github.com/genid/MatchY/issues)
2. Click "New Issue"
3. Select "Bug Report"
4. Fill out all sections:
   - Clear description of the bug
   - Steps to reproduce
   - Expected behavior
   - Actual behavior
   - Environment details (OS, Python version, etc.)
   - Screenshots if applicable

**Before reporting:**
- Check if the issue already exists
- Verify the bug is reproducible
- Check if it's fixed in the latest version

## Suggesting Features

Have an idea for a new feature?

1. Go to [Issues](https://github.com/genid/MatchY/issues)
2. Click "New Issue"
3. Select "Feature Request"
4. Describe:
   - The problem or use case
   - Proposed solution
   - Alternative solutions considered
   - Why this would be useful

## Questions?

- Create a [Question Issue](https://github.com/genid/MatchY/issues)
- Check existing issues for similar questions
- Read the [User Manual](USER_MANUAL.md)

## License

By contributing to MatchY, you agree that your contributions will be licensed under the MIT License.

## Acknowledgments

Thank you for contributing to MatchY! Your efforts help make this tool better for the forensic genetics community.
