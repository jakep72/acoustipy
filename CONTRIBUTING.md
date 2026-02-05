# Contributing to acoustipy

Thank you for your interest in contributing to acoustipy! This document provides guidelines and instructions for contributing.

## Development Setup

### Prerequisites

- Python 3.9 or higher
- Git

### Setting Up the Development Environment

1. **Clone the repository:**
   ```bash
   git clone https://github.com/jakep72/acoustipy.git
   cd acoustipy
   ```

2. **Create and activate a virtual environment:**
   ```bash
   python -m venv .venv
   source .venv/bin/activate  # Linux/macOS
   # or
   .venv\Scripts\activate  # Windows
   ```

3. **Install the package in editable mode with dev dependencies:**
   ```bash
   pip install -e ".[dev]"
   ```

   Or with uv (recommended for faster installation):
   ```bash
   uv pip install -e ".[dev]"
   ```

### GPU Support (Optional)

The default installation uses CPU-only PyTorch. For GPU acceleration during development:

```bash
# Install PyTorch with CUDA support first
pip install torch --index-url https://download.pytorch.org/whl/cu121

# Then install acoustipy
pip install -e ".[dev]"
```

See [PyTorch's installation guide](https://pytorch.org/get-started/locally/) for your specific CUDA version.

## Running Tests

### Run all tests:
```bash
pytest tests/ -v
```

### Run with coverage:
```bash
pytest tests/ --cov=acoustipy --cov-report=html
```

### Run a specific test file:
```bash
pytest tests/test_tmm.py -v
```

### Run a specific test:
```bash
pytest tests/test_paramfinder.py::test_inverse -v
```

## Code Style

- Follow PEP 8 style guidelines
- Use type hints for function parameters and return values
- Add docstrings in NumPy format for all public methods
- Keep lines under 100 characters where practical

### Example docstring format:
```python
def my_function(param1: float, param2: str = "default") -> dict:
    """
    Short description of the function.

    Longer description if needed.

    Parameters
    ----------
    param1 : float
        Description of param1.
    param2 : str, optional
        Description of param2 (default is "default").

    Returns
    -------
    result : dict
        Description of the return value.

    Raises
    ------
    ValueError
        If param1 is negative.

    Examples
    --------
    >>> my_function(1.5, "test")
    {'key': 'value'}
    """
    pass
```

## Making Changes

1. **Create a new branch:**
   ```bash
   git checkout -b feature/your-feature-name
   ```

2. **Make your changes** and ensure tests pass:
   ```bash
   pytest tests/ -v
   ```

3. **Commit your changes:**
   ```bash
   git add .
   git commit -m "Add feature: description of changes"
   ```

4. **Push to your fork and create a pull request:**
   ```bash
   git push origin feature/your-feature-name
   ```

## Pull Request Guidelines

- Provide a clear description of the changes
- Include tests for new functionality
- Update documentation if needed
- Ensure all tests pass
- Keep commits focused and atomic

## Project Structure

```
acoustipy/
├── src/
│   └── acoustipy/
│       ├── __init__.py      # Package exports
│       ├── TMM.py           # AcousticTMM class (transfer matrix method)
│       ├── Params.py        # AcousticID class (parameter identification)
│       └── Database.py      # AcoustiBase class (database operations)
├── tests/
│   ├── test_tmm.py          # Tests for AcousticTMM
│   └── test_paramfinder.py  # Tests for AcousticID
├── docs/                    # Documentation source
├── pyproject.toml           # Package configuration
└── README.md
```

## Key Classes

### AcousticTMM
The acoustic transfer matrix method implementation. Handles:
- Layer definitions (JCA, JCAL, JCAPL, Biot, MPP, etc.)
- Transfer matrix calculations
- Absorption, reflection, and transmission coefficients

### AcousticID
Parameter identification routines:
- **Inverse**: Gradient-based optimization using SLSQP
- **Indirect**: Analytical parameter estimation from impedance data
- **Hybrid**: Combines Inverse and Indirect methods
- **ML**: Machine learning approach using Adam optimizer

## Reporting Issues

When reporting issues, please include:
- Python version
- Operating system
- Steps to reproduce the issue
- Expected vs. actual behavior
- Any error messages or tracebacks

## Questions?

Feel free to open an issue for any questions about contributing.
