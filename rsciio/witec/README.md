# WiTec WIP File Format Support

This plugin provides support for reading WiTec WIP (WiTec Project) files. WiTec WIP files are binary files used by WiTec's Project software to store confocal Raman microscopy data.

## Supported Features

- **Spectral Maps**: 3D datasets containing Raman spectra for each spatial position
- **Images**: 2D grayscale or color images from optical or other microscopy modes
- **Individual Spectra**: Single Raman spectra
- **Spatial Calibration**: Pixel size and positioning information
- **Spectral Calibration**: Wavenumber axis calibration
- **Fallback Raw Data Extraction**: For files that don't match known structures

## File Format

WiTec WIP files typically contain:

1. **File Header**: Contains format version and dataset information
2. **Data Blocks**: Structured blocks containing different data types:
   - `DATA`: Spectral map data with 3D arrays
   - `IMG `: Image data with 2D arrays  
   - `SPEC`: Individual spectrum data
   - `META`: Metadata information
   - `AXIS`: Axis calibration data
   - `COOR`: Coordinate information

## Usage

```python
import rsciio

# Read WiTec WIP file
data = rsciio.load("spectrum.wip")

# Access the data
spectra = data[0].data  # 3D array for spectral maps: (y, x, spectral_points)
axes = data[0].axes     # Axis information with calibration

# For spectral maps:
# axes[0]: y-axis (spatial, µm)
# axes[1]: x-axis (spatial, µm) 
# axes[2]: Raman shift axis (cm⁻¹)
```

## Metadata

The plugin preserves:
- **Spatial calibration**: Pixel sizes and offsets
- **Spectral calibration**: Wavenumber range and scaling
- **Original metadata**: Raw file header and block information
- **Signal type**: Automatically set to 'Raman' for spectral data or 'Image' for images

## Limitations

- Lazy loading is not currently supported
- Some proprietary WiTec metadata may not be fully parsed
- File format reverse-engineered from typical scientific data patterns

## Notes

This implementation is based on analysis of WiTec file structures and common patterns in scientific data formats. The actual WiTec WIP format may vary between software versions, so some files might fall back to raw binary extraction if structured parsing fails.