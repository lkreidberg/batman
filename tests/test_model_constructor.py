import numpy as np
import pytest
from batman import TransitModel, TransitParams  # Adjust the import according to your project structure


def test_valid_initialization():
    """Test the constructor with valid input parameters."""
    params = TransitParams()
    params.limb_dark = "quadratic"
    params.t0 = 0.0
    params.per = 1.0
    params.rp = 0.1
    params.a = 15.0
    params.inc = 87.0
    params.ecc = 0.0
    params.w = 90.0
    params.u = [0.1, 0.3]
    t = np.linspace(-0.015, 0.015, 100)

    model = TransitModel(params, t)
    assert model is not None
    assert model.t.size == 100


def test_invalid_coefficients():
    """Test the constructor with an invalid number of limb darkening coefficients."""
    params = TransitParams()
    params.limb_dark = "quadratic"
    params.u = [0.1]  # Incorrect number for quadratic
    params.rp = 0.1

    t = np.linspace(0, 10, 100)
    with pytest.raises(Exception) as excinfo:
        TransitModel(params, t)
    assert "Incorrect number of coefficients" in str(excinfo.value)


def test_unsupported_limb_darkening():
    """Test with an unsupported limb darkening model."""
    params = TransitParams()
    params.limb_dark = "unsupported_model"
    params.u = [0.1, 0.2]
    params.rp = 0.1

    t = np.linspace(0, 10, 100)
    with pytest.raises(Exception) as excinfo:
        TransitModel(params, t)
    assert "limb darkening not supported" in str(excinfo.value)


def test_invalid_error_tolerance():
    """Test with an invalid error tolerance."""
    params = TransitParams()
    params.limb_dark = "quadratic"
    params.u = [0.1, 0.2]
    params.rp = 0.1

    t = np.linspace(0, 10, 100)
    with pytest.raises(Exception) as excinfo:
        TransitModel(params, t, max_err=0.0001)  # Below the threshold
    assert "The lowest allowed value for max_err is 0.001" in str(excinfo.value)


def test_invalid_transit_type():
    """Test with an invalid transit type."""
    params = TransitParams()
    params.limb_dark = "quadratic"
    params.u = [0.1, 0.2]
    params.rp = 0.1

    t = np.linspace(0, 10, 100)
    with pytest.raises(Exception) as excinfo:
        TransitModel(params, t, transittype="invalid_type")
    assert "Allowed transit types are \"primary\" and \"secondary\"" in str(excinfo.value)


def test_invalid_exposure_time_with_supersampling():
    """Test with an invalid exposure time when supersampling is used."""
    params = TransitParams()
    params.limb_dark = "quadratic"
    params.u = [0.1, 0.2]
    params.rp = 0.1

    t = np.linspace(0, 10, 100)
    with pytest.raises(Exception) as excinfo:
        TransitModel(params, t, supersample_factor=2,
                     exp_time=0)  # Invalid because exp_time must be > 0 with supersampling
    assert "Please enter a valid exposure time" in str(excinfo.value)


def test_non_array_time_input():
    """Test with a non-array time input."""
    params = TransitParams()
    params.limb_dark = "quadratic"
    params.u = [0.1, 0.2]
    params.rp = 0.1

    t = list(np.linspace(0, 10, 100))  # Pass a list instead of a numpy array
    with pytest.raises(Exception) as excinfo:
        TransitModel(params, t)
    assert "Times t must be a numpy array" in str(excinfo.value)

# Use the following command in your terminal to run these tests:
# pytest test_transit_model_constructor.py
