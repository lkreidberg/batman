import numpy as np
import pytest
from batman import TransitModel


def test_valid_initialization_quadratic(default_params, default_time):
    """Test the constructor with valid input parameters."""
    model = TransitModel(default_params, default_time)
    assert model is not None
    assert model.t.size == 100


def test_valid_initialization_linear(default_params, default_time):
    """Test the constructor with valid input parameters."""
    default_params.limb_dark = "linear"
    default_params.u = [0.1]
    model = TransitModel(default_params, default_time)
    assert model is not None
    assert model.t.size == 100


def test_valid_initialization_nonlinear(default_params, default_time):
    """Test the constructor with valid input parameters."""
    default_params.limb_dark = "nonlinear"
    default_params.u = [0.1, 0.3, 0.2, 0.1]
    model = TransitModel(default_params, default_time)
    assert model is not None
    assert model.t.size == 100


def test_valid_initialization_squareroot(default_params, default_time):
    """Test the constructor with valid input parameters."""
    default_params.limb_dark = "squareroot"
    model = TransitModel(default_params, default_time)
    assert model is not None
    assert model.t.size == 100


def test_valid_initialization_uniform(default_params, default_time):
    """Test the constructor with valid input parameters."""
    default_params.limb_dark = "uniform"
    default_params.u = []
    model = TransitModel(default_params, default_time)
    assert model is not None
    assert model.t.size == 100


def test_valid_initialization_logarithmic(default_params, default_time):
    """Test the constructor with valid input parameters."""
    default_params.limb_dark = "logarithmic"
    model = TransitModel(default_params, default_time)
    assert model is not None
    assert model.t.size == 100


def test_valid_initialization_exponential(default_params, default_time):
    """Test the constructor with valid input parameters."""
    default_params.limb_dark = "exponential"
    model = TransitModel(default_params, default_time)
    assert model is not None
    assert model.t.size == 100


def test_valid_initialization_power2(default_params, default_time):
    """Test the constructor with valid input parameters."""
    default_params.limb_dark = "power2"
    model = TransitModel(default_params, default_time)
    assert model is not None
    assert model.t.size == 100


def test_valid_initialization_custom(default_params, default_time):
    """Test the constructor with valid input parameters."""
    default_params.limb_dark = "custom"
    default_params.u = [0.1, 0.3, 0.2, 0.1, 0.2, 0.3]
    model = TransitModel(default_params, default_time)
    assert model is not None
    assert model.t.size == 100


def test_invalid_coefficients(default_params, default_time):
    """Test the constructor with an invalid number of limb darkening coefficients."""
    default_params.u = [0.1]
    with pytest.raises(Exception) as excinfo:
        TransitModel(default_params, default_time)
    assert "Incorrect number of coefficients" in str(excinfo.value)


def test_unsupported_limb_darkening(default_params, default_time):
    """Test with an unsupported limb darkening model."""
    default_params.limb_dark = "unsupported_model"
    with pytest.raises(Exception) as excinfo:
        TransitModel(default_params, default_time)
    assert "limb darkening not supported" in str(excinfo.value)


def test_invalid_error_tolerance(default_params, default_time):
    """Test with an invalid error tolerance."""
    with pytest.raises(Exception) as excinfo:
        TransitModel(
            default_params, default_time, max_err=0.0001
        )  # Below the threshold
    assert "The lowest allowed value for max_err is 0.001" in str(excinfo.value)


def test_invalid_transit_type(default_params, default_time):
    """Test with an invalid transit type."""
    with pytest.raises(Exception) as excinfo:
        TransitModel(default_params, default_time, transittype="invalid_type")
    assert 'Allowed transit types are "primary" and "secondary"' in str(excinfo.value)


def test_invalid_exposure_time_with_supersampling(default_params, default_time):
    """Test with an invalid exposure time when supersampling is used."""
    with pytest.raises(Exception) as excinfo:
        TransitModel(
            default_params, default_time, supersample_factor=2, exp_time=0
        )  # Invalid because exp_time must be > 0 with supersampling
    assert "Please enter a valid exposure time" in str(excinfo.value)


def test_non_array_time_input(default_params):
    t = list(np.linspace(0, 10, 100))  # Pass a list instead of a numpy array
    with pytest.raises(Exception) as excinfo:
        TransitModel(default_params, t)
    assert "Times t must be a numpy array" in str(excinfo.value)
