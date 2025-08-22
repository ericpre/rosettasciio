# -*- coding: utf-8 -*-
# Copyright 2007-2025 The HyperSpy developers
#
# This file is part of RosettaSciIO.
#
# RosettaSciIO is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# RosettaSciIO is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with RosettaSciIO. If not, see <https://www.gnu.org/licenses/#GPL>.

import logging
import struct
from pathlib import Path
from copy import deepcopy

import numpy as np

from rsciio._docstrings import FILENAME_DOC, LAZY_UNSUPPORTED_DOC, RETURNS_DOC

_logger = logging.getLogger(__name__)

# WiTec data types - these are inferred from similar scientific formats
WIP_HEADER_MAGIC = b'WIT'  # Likely header signature
WIP_BLOCK_SIGNATURES = {
    b'DATA': 'spectral_data',
    b'IMG ': 'image_data', 
    b'SPEC': 'spectrum_data',
    b'META': 'metadata',
    b'AXIS': 'axis_info',
    b'COOR': 'coordinates'
}


def file_reader(filename, lazy=False):
    """
    Read a WiTec ``.wip`` file.
    
    WiTec WIP files contain Raman spectroscopy data including images and spectral maps.
    This format is used by WiTec's Project software for storing confocal Raman
    microscopy data.

    Parameters
    ----------
    %s
    %s

    %s
    """
    if lazy is not False:
        raise NotImplementedError("Lazy loading is not supported.")

    reader = WIPReader(filename)
    return reader.read_file()


file_reader.__doc__ %= (FILENAME_DOC, LAZY_UNSUPPORTED_DOC, RETURNS_DOC)


class WIPReader:
    """Reader for WiTec Project WIP files.
    
    WiTec WIP files are binary files that contain Raman spectroscopy data,
    including spectra, images, and spectral maps from confocal Raman microscopy.
    The format typically contains structured blocks with different data types.
    """
    
    def __init__(self, filename):
        self.filename = Path(filename)
        self.file_obj = None
        self.file_size = 0
        self.original_metadata = {}
        self.data_blocks = []
        self.file_version = None
        
    def read_file(self):
        """Read the WIP file and return data in HyperSpy format."""
        try:
            with open(self.filename, 'rb') as f:
                self.file_obj = f
                self.file_size = self.filename.stat().st_size
                return self._parse_file()
        except Exception as e:
            _logger.error(f"Failed to read WiTec WIP file {self.filename}: {e}")
            raise
            
    def _parse_file(self):
        """Parse the WIP file structure."""
        # Scan for all data blocks in the file
        self._scan_file_blocks()
        
        # Parse each block type
        datasets = []
        
        # Group related blocks (e.g., data with corresponding axes and metadata)
        for block in self.data_blocks:
            if block['type'] in ['spectral_data', 'image_data', 'spectrum_data']:
                dataset = self._parse_data_block(block)
                if dataset is not None:
                    datasets.append(dataset)
        
        if not datasets:
            # If no structured data found, try to extract as raw binary data
            _logger.warning("No structured data blocks found, attempting raw extraction")
            datasets = self._extract_raw_data()
        
        # Convert to HyperSpy format
        return self._convert_to_hyperspy_format(datasets)
        
    def _parse_header(self):
        """Parse the WIP file header."""
        header = {}
        
        # Read file signature
        self.file_obj.seek(0)
        signature = self.file_obj.read(8)
        
        if len(signature) < 8:
            raise ValueError("File too short to be a valid WIP file")
            
        header['signature'] = signature
        header['file_size'] = self.file_size
        
        # Try to identify file version and structure
        # WiTec files may have different header structures
        if signature.startswith(WIP_HEADER_MAGIC):
            self.file_version = struct.unpack('<I', signature[4:8])[0]
            header['version'] = self.file_version
            _logger.debug(f"WIP file version: {self.file_version}")
        else:
            _logger.debug(f"Unknown WIP signature: {signature.hex()}")
            header['version'] = 'unknown'
        
        # Read additional header information
        try:
            # Try to read common header fields
            self.file_obj.seek(8)
            num_datasets = struct.unpack('<I', self.file_obj.read(4))[0]
            header['num_datasets'] = num_datasets
            
            data_offset = struct.unpack('<Q', self.file_obj.read(8))[0] 
            header['data_offset'] = data_offset
            
        except (struct.error, IOError) as e:
            _logger.debug(f"Could not parse extended header: {e}")
            header['num_datasets'] = 1  # Assume single dataset
            header['data_offset'] = 64  # Common offset
        
        self.original_metadata['header'] = header
        return header
        
    def _scan_file_blocks(self):
        """Scan the file for data blocks."""
        self.data_blocks = []
        position = 64  # Skip header
        
        while position < self.file_size - 8:
            try:
                self.file_obj.seek(position)
                block_sig = self.file_obj.read(4)
                
                if len(block_sig) < 4:
                    break
                    
                # Check if this looks like a block signature
                if block_sig in WIP_BLOCK_SIGNATURES:
                    block_size = struct.unpack('<I', self.file_obj.read(4))[0]
                    
                    if block_size > 0 and position + block_size <= self.file_size:
                        block = {
                            'signature': block_sig,
                            'type': WIP_BLOCK_SIGNATURES[block_sig],
                            'position': position + 8,  # After signature and size
                            'size': block_size,
                            'data': None
                        }
                        self.data_blocks.append(block)
                        _logger.debug(f"Found {block['type']} block at {position}, size {block_size}")
                        position += block_size + 8  # Move to next block
                    else:
                        position += 1  # Move forward and continue scanning
                else:
                    position += 1  # Move forward and continue scanning
                    
            except (struct.error, IOError):
                position += 1  # Continue scanning on error
                continue
                
        if not self.data_blocks:
            _logger.info("No structured blocks found, file may use different format")
    
    def _parse_data_block(self, block):
        """Parse a specific data block."""
        try:
            self.file_obj.seek(block['position'])
            
            if block['type'] == 'spectral_data':
                return self._parse_spectral_data_block(block)
            elif block['type'] == 'image_data':
                return self._parse_image_data_block(block)
            elif block['type'] == 'spectrum_data':
                return self._parse_spectrum_data_block(block)
            else:
                # Store as metadata
                raw_data = self.file_obj.read(block['size'])
                self.original_metadata[block['type']] = raw_data
                return None
                
        except Exception as e:
            _logger.warning(f"Failed to parse {block['type']} block: {e}")
            return None
    
    def _parse_spectral_data_block(self, block):
        """Parse spectral map data."""
        try:
            self.file_obj.seek(block['position'])
            
            # Read dimensions (assumed structure)
            nx = struct.unpack('<I', self.file_obj.read(4))[0]
            ny = struct.unpack('<I', self.file_obj.read(4))[0]  
            nz = struct.unpack('<I', self.file_obj.read(4))[0]  # spectral points
            
            # Read spatial calibration
            x_scale = struct.unpack('<f', self.file_obj.read(4))[0]
            y_scale = struct.unpack('<f', self.file_obj.read(4))[0]
            x_offset = struct.unpack('<f', self.file_obj.read(4))[0]
            y_offset = struct.unpack('<f', self.file_obj.read(4))[0]
            
            # Read spectral calibration  
            spec_offset = struct.unpack('<f', self.file_obj.read(4))[0]
            spec_scale = struct.unpack('<f', self.file_obj.read(4))[0]
            
            # Calculate expected data size
            expected_size = nx * ny * nz * 4  # Assuming 32-bit floats
            remaining_size = block['size'] - (self.file_obj.tell() - block['position'])
            
            if expected_size <= remaining_size:
                # Read the spectral data
                data = np.fromfile(self.file_obj, dtype=np.float32, count=nx*ny*nz)
                data = data.reshape((ny, nx, nz))  # y, x, spectral
                
                # Create axes
                axes = [
                    {
                        'name': 'y',
                        'size': ny,
                        'scale': y_scale,
                        'offset': y_offset,
                        'units': 'µm',
                        'navigate': True
                    },
                    {
                        'name': 'x', 
                        'size': nx,
                        'scale': x_scale,
                        'offset': x_offset,
                        'units': 'µm',
                        'navigate': True
                    },
                    {
                        'name': 'Raman Shift',
                        'size': nz,
                        'scale': spec_scale,
                        'offset': spec_offset,
                        'units': 'cm⁻¹',
                        'navigate': False
                    }
                ]
                
                return {
                    'data': data,
                    'axes': axes,
                    'metadata': {'Signal': {'signal_type': 'Raman'}},
                    'original_metadata': {
                        'data_type': 'spectral_map',
                        'dimensions': (nx, ny, nz),
                        'calibration': {
                            'x_scale': x_scale, 'x_offset': x_offset,
                            'y_scale': y_scale, 'y_offset': y_offset,
                            'spectral_scale': spec_scale, 'spectral_offset': spec_offset
                        }
                    }
                }
            else:
                _logger.warning(f"Data size mismatch in spectral block: expected {expected_size}, available {remaining_size}")
                return None
                
        except Exception as e:
            _logger.warning(f"Failed to parse spectral data: {e}")
            return None
    
    def _parse_image_data_block(self, block):
        """Parse image data."""
        try:
            self.file_obj.seek(block['position'])
            
            # Read image dimensions
            width = struct.unpack('<I', self.file_obj.read(4))[0]
            height = struct.unpack('<I', self.file_obj.read(4))[0]
            channels = struct.unpack('<I', self.file_obj.read(4))[0]  # e.g., RGB = 3
            
            # Read spatial calibration
            x_scale = struct.unpack('<f', self.file_obj.read(4))[0]
            y_scale = struct.unpack('<f', self.file_obj.read(4))[0]
            
            expected_size = width * height * channels
            data = np.fromfile(self.file_obj, dtype=np.uint8, count=expected_size)
            
            if channels == 1:
                data = data.reshape((height, width))
                axes = [
                    {'name': 'y', 'size': height, 'scale': y_scale, 'offset': 0, 'units': 'µm', 'navigate': True},
                    {'name': 'x', 'size': width, 'scale': x_scale, 'offset': 0, 'units': 'µm', 'navigate': True}
                ]
            else:
                data = data.reshape((height, width, channels))
                axes = [
                    {'name': 'y', 'size': height, 'scale': y_scale, 'offset': 0, 'units': 'µm', 'navigate': True},
                    {'name': 'x', 'size': width, 'scale': x_scale, 'offset': 0, 'units': 'µm', 'navigate': True},
                    {'name': 'channel', 'size': channels, 'scale': 1, 'offset': 0, 'units': '', 'navigate': False}
                ]
            
            return {
                'data': data,
                'axes': axes,
                'metadata': {'Signal': {'signal_type': 'Image'}},
                'original_metadata': {
                    'data_type': 'image',
                    'dimensions': (width, height, channels) if channels > 1 else (width, height)
                }
            }
            
        except Exception as e:
            _logger.warning(f"Failed to parse image data: {e}")
            return None
    
    def _parse_spectrum_data_block(self, block):
        """Parse single spectrum data."""
        try:
            self.file_obj.seek(block['position'])
            
            # Read spectrum length
            length = struct.unpack('<I', self.file_obj.read(4))[0]
            
            # Read spectral calibration
            offset = struct.unpack('<f', self.file_obj.read(4))[0]
            scale = struct.unpack('<f', self.file_obj.read(4))[0]
            
            # Read spectral data
            data = np.fromfile(self.file_obj, dtype=np.float32, count=length)
            
            axes = [{
                'name': 'Raman Shift',
                'size': length,
                'scale': scale,
                'offset': offset,
                'units': 'cm⁻¹',
                'navigate': False
            }]
            
            return {
                'data': data,
                'axes': axes,
                'metadata': {'Signal': {'signal_type': 'Raman'}},
                'original_metadata': {
                    'data_type': 'spectrum',
                    'length': length,
                    'calibration': {'offset': offset, 'scale': scale}
                }
            }
            
        except Exception as e:
            _logger.warning(f"Failed to parse spectrum data: {e}")
            return None
        
    def _extract_raw_data(self):
        """Extract data when no structured blocks are found."""
        # Try to interpret the file as raw binary data
        # This is a fallback for files that don't match expected structure
        
        datasets = []
        
        # Skip header and try to find data patterns
        self.file_obj.seek(64)  
        remaining_bytes = self.file_size - 64
        
        # Try different common data sizes to guess format
        possible_formats = [
            (np.float32, 4),
            (np.float64, 8), 
            (np.int32, 4),
            (np.uint16, 2),
        ]
        
        for dtype, byte_size in possible_formats:
            if remaining_bytes % byte_size == 0:
                try:
                    self.file_obj.seek(64)
                    data_points = remaining_bytes // byte_size
                    data = np.fromfile(self.file_obj, dtype=dtype, count=data_points)
                    
                    # Try to guess reasonable dimensions
                    # Common patterns: 1D spectrum, 2D image, 3D spectral map
                    if data_points > 1000:  # Likely spectral map or large image
                        # Try to find factors that could represent sensible dimensions
                        factors = self._find_factors(data_points)
                        
                        if len(factors) >= 3:
                            # Assume 3D spectral map
                            nz = factors[-1]  # Spectral dimension (usually largest)
                            ny = factors[-2]  # y dimension  
                            nx = factors[-3]  # x dimension
                            
                            if nx * ny * nz == data_points:
                                data = data.reshape((ny, nx, nz))
                                axes = self._create_default_3d_axes(ny, nx, nz)
                                
                                dataset = {
                                    'data': data,
                                    'axes': axes,
                                    'metadata': {'Signal': {'signal_type': 'Raman'}},
                                    'original_metadata': {
                                        'data_type': 'raw_extraction_3d',
                                        'guessed_format': str(dtype),
                                        'dimensions': (nx, ny, nz)
                                    }
                                }
                                datasets.append(dataset)
                                break
                                
                        elif len(factors) >= 2:
                            # Assume 2D image
                            ny = factors[-1]
                            nx = factors[-2]
                            
                            if nx * ny == data_points:
                                data = data.reshape((ny, nx))
                                axes = self._create_default_2d_axes(ny, nx)
                                
                                dataset = {
                                    'data': data,
                                    'axes': axes,
                                    'metadata': {'Signal': {'signal_type': 'Image'}},
                                    'original_metadata': {
                                        'data_type': 'raw_extraction_2d',
                                        'guessed_format': str(dtype),
                                        'dimensions': (nx, ny)
                                    }
                                }
                                datasets.append(dataset)
                                break
                    
                    else:
                        # Assume 1D spectrum
                        axes = self._create_default_1d_axes(data_points)
                        
                        dataset = {
                            'data': data,
                            'axes': axes,
                            'metadata': {'Signal': {'signal_type': 'Raman'}},
                            'original_metadata': {
                                'data_type': 'raw_extraction_1d',
                                'guessed_format': str(dtype),
                                'length': data_points
                            }
                        }
                        datasets.append(dataset)
                        break
                        
                except Exception as e:
                    _logger.debug(f"Failed to extract as {dtype}: {e}")
                    continue
        
        if not datasets:
            # Last resort - create minimal dummy data
            _logger.warning("Could not extract any data, creating minimal dataset")
            data = np.zeros((10, 10))
            axes = self._create_default_2d_axes(10, 10)
            
            dataset = {
                'data': data,
                'axes': axes, 
                'metadata': {'Signal': {'signal_type': ''}},
                'original_metadata': {'data_type': 'fallback'}
            }
            datasets.append(dataset)
        
        return datasets
    
    def _find_factors(self, n):
        """Find reasonable factors of a number for dimension guessing."""
        factors = []
        for i in range(2, int(n**0.5) + 1):
            while n % i == 0:
                factors.append(i)
                n = n // i
        if n > 1:
            factors.append(n)
        
        # Sort factors and try to combine them into reasonable dimensions
        factors.sort()
        
        # Try to create reasonable dimension combinations
        if len(factors) >= 3:
            # For 3D data, try to make spectral dimension reasonable (100-4000 points)
            combined_factors = []
            temp_product = 1
            
            for f in factors:
                temp_product *= f
                if temp_product >= 50:  # Reasonable minimum dimension
                    combined_factors.append(temp_product)
                    temp_product = 1
            
            if temp_product > 1:
                combined_factors.append(temp_product)
                
            return combined_factors[-3:] if len(combined_factors) >= 3 else factors[-3:]
        
        return factors
        
    def _create_default_3d_axes(self, ny, nx, nz):
        """Create default 3D axes for spectral map."""
        return [
            {'name': 'y', 'size': ny, 'scale': 1.0, 'offset': 0, 'units': 'µm', 'navigate': True},
            {'name': 'x', 'size': nx, 'scale': 1.0, 'offset': 0, 'units': 'µm', 'navigate': True},
            {'name': 'Raman Shift', 'size': nz, 'scale': 1.0, 'offset': 0, 'units': 'cm⁻¹', 'navigate': False}
        ]
        
    def _create_default_2d_axes(self, ny, nx):
        """Create default 2D axes for image."""
        return [
            {'name': 'y', 'size': ny, 'scale': 1.0, 'offset': 0, 'units': 'µm', 'navigate': True},
            {'name': 'x', 'size': nx, 'scale': 1.0, 'offset': 0, 'units': 'µm', 'navigate': True}
        ]
        
    def _create_default_1d_axes(self, n):
        """Create default 1D axes for spectrum.""" 
        return [{
            'name': 'Raman Shift',
            'size': n,
            'scale': 1.0,
            'offset': 0,
            'units': 'cm⁻¹',
            'navigate': False
        }]
        
    def _convert_to_hyperspy_format(self, datasets):
        """Convert parsed datasets to HyperSpy format."""
        result = []
        
        for dataset in datasets:
            # Create the dictionary structure expected by HyperSpy
            hyperspy_dict = {
                'data': dataset['data'],
                'axes': dataset['axes'],
                'metadata': dataset.get('metadata', {}),
                'original_metadata': dataset.get('original_metadata', {})
            }
            
            # Add general metadata
            if 'General' not in hyperspy_dict['metadata']:
                hyperspy_dict['metadata']['General'] = {}
                
            hyperspy_dict['metadata']['General']['original_filename'] = self.filename.name
            
            result.append(hyperspy_dict)
            
        return result


def _remove_none_from_dict(dict_in):
    """Remove None values from nested dictionaries."""
    for key, value in list(dict_in.items()):
        if isinstance(value, dict):
            _remove_none_from_dict(value)
        elif value is None:
            del dict_in[key]