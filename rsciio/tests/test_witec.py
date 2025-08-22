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

import os
import struct
from pathlib import Path

import numpy as np
import pytest

from rsciio.witec import file_reader


class TestWiTecReader:
    @pytest.fixture
    def spectral_map_file(self, tmp_path):
        """Create a test spectral map file."""
        filename = tmp_path / "test_spectral.wip"
        
        with open(filename, 'wb') as f:
            # Write file header
            f.write(b'WIT\x00')  # Magic signature
            f.write(struct.pack('<I', 1))  # Version 1
            f.write(struct.pack('<I', 1))  # Number of datasets 
            f.write(struct.pack('<Q', 64))  # Data offset
            f.write(b'\x00' * 44)  # Padding to 64 bytes
            
            # Write a spectral data block
            f.write(b'DATA')  # Block signature
            nx, ny, nz = 5, 4, 50
            block_size = 4*3 + 4*6 + 4*nx*ny*nz  # dimensions + calibration + data
            f.write(struct.pack('<I', block_size))
            
            # Dimensions
            f.write(struct.pack('<I', nx))  # nx
            f.write(struct.pack('<I', ny))  # ny  
            f.write(struct.pack('<I', nz))  # nz spectral
            
            # Spatial calibration
            f.write(struct.pack('<f', 0.5))  # x_scale
            f.write(struct.pack('<f', 0.4))  # y_scale
            f.write(struct.pack('<f', 1.0))  # x_offset
            f.write(struct.pack('<f', 2.0))  # y_offset
            
            # Spectral calibration  
            f.write(struct.pack('<f', 100.0))  # spec_offset
            f.write(struct.pack('<f', 5.0))    # spec_scale
            
            # Generate test data
            data = np.random.random((ny, nx, nz)).astype(np.float32) * 1000
            data.tofile(f)
        
        return filename

    @pytest.fixture
    def image_file(self, tmp_path):
        """Create a test image file."""
        filename = tmp_path / "test_image.wip"
        
        with open(filename, 'wb') as f:
            # Write file header
            f.write(b'WIT\x00')  # Magic signature
            f.write(struct.pack('<I', 1))  # Version 1
            f.write(struct.pack('<I', 1))  # Number of datasets 
            f.write(struct.pack('<Q', 64))  # Data offset
            f.write(b'\x00' * 44)  # Padding to 64 bytes
            
            # Write image data block
            f.write(b'IMG ')  # Block signature
            width, height, channels = 32, 24, 1
            block_size = 4*3 + 4*2 + width*height*channels
            f.write(struct.pack('<I', block_size))
            
            # Image dimensions
            f.write(struct.pack('<I', width))
            f.write(struct.pack('<I', height)) 
            f.write(struct.pack('<I', channels))
            
            # Spatial calibration
            f.write(struct.pack('<f', 0.1))  # x_scale
            f.write(struct.pack('<f', 0.1))  # y_scale
            
            # Generate test image
            image_data = (np.random.random((height, width)) * 255).astype(np.uint8)
            image_data.tofile(f)
        
        return filename

    @pytest.fixture
    def raw_file(self, tmp_path):
        """Create a raw binary file for fallback testing."""
        filename = tmp_path / "test_raw.wip"
        
        with open(filename, 'wb') as f:
            # Write minimal header
            f.write(b'UNKNOWN\x00')
            f.write(b'\x00' * 56)  # Padding to 64 bytes
            
            # Write some float data that should be interpreted as 1D spectrum
            data = np.random.random(200).astype(np.float32) * 1000
            data.tofile(f)
        
        return filename

    def test_spectral_map_reading(self, spectral_map_file):
        """Test reading spectral map data."""
        results = file_reader(spectral_map_file)
        
        assert len(results) == 1
        dataset = results[0]
        
        # Check data shape and type
        assert dataset['data'].shape == (4, 5, 50)  # ny, nx, nz
        assert dataset['original_metadata']['data_type'] == 'spectral_map'
        assert dataset['metadata']['Signal']['signal_type'] == 'Raman'
        
        # Check axes
        assert len(dataset['axes']) == 3
        assert dataset['axes'][0]['name'] == 'y'
        assert dataset['axes'][1]['name'] == 'x'
        assert dataset['axes'][2]['name'] == 'Raman Shift'
        
        # Check calibration
        assert dataset['axes'][1]['scale'] == pytest.approx(0.5, abs=1e-6)
        assert dataset['axes'][0]['scale'] == pytest.approx(0.4, abs=1e-6)
        assert dataset['axes'][2]['scale'] == pytest.approx(5.0, abs=1e-6)
        assert dataset['axes'][2]['offset'] == pytest.approx(100.0, abs=1e-6)

    def test_image_reading(self, image_file):
        """Test reading image data."""
        results = file_reader(image_file)
        
        assert len(results) == 1
        dataset = results[0]
        
        # Check data shape and type
        assert dataset['data'].shape == (24, 32)  # height, width
        assert dataset['original_metadata']['data_type'] == 'image'
        assert dataset['metadata']['Signal']['signal_type'] == 'Image'
        
        # Check axes
        assert len(dataset['axes']) == 2
        assert dataset['axes'][0]['name'] == 'y'
        assert dataset['axes'][1]['name'] == 'x'
        
        # Check calibration
        assert dataset['axes'][0]['scale'] == pytest.approx(0.1, abs=1e-6)
        assert dataset['axes'][1]['scale'] == pytest.approx(0.1, abs=1e-6)

    def test_raw_data_fallback(self, raw_file):
        """Test fallback raw data extraction."""
        results = file_reader(raw_file)
        
        assert len(results) == 1
        dataset = results[0]
        
        # Check that data was extracted
        assert dataset['data'].size == 200
        assert 'raw_extraction' in dataset['original_metadata']['data_type']
        
        # Should default to 1D spectrum for small data
        assert dataset['metadata']['Signal']['signal_type'] == 'Raman'
        assert len(dataset['axes']) == 1
        assert dataset['axes'][0]['name'] == 'Raman Shift'

    def test_invalid_file(self, tmp_path):
        """Test handling of invalid files.""" 
        filename = tmp_path / "empty.wip"
        
        # Create empty file
        with open(filename, 'wb') as f:
            pass
        
        with pytest.raises(ValueError, match="File too short"):
            file_reader(filename)

    def test_lazy_loading_not_supported(self, spectral_map_file):
        """Test that lazy loading raises NotImplementedError."""
        with pytest.raises(NotImplementedError, match="Lazy loading is not supported"):
            file_reader(spectral_map_file, lazy=True)

    def test_metadata_preservation(self, spectral_map_file):
        """Test that metadata is properly preserved."""
        results = file_reader(spectral_map_file)
        dataset = results[0]
        
        # Check that original filename is preserved
        assert 'General' in dataset['metadata']
        assert 'original_filename' in dataset['metadata']['General']
        assert dataset['metadata']['General']['original_filename'] == 'test_spectral.wip'
        
        # Check original metadata structure
        assert 'header' in dataset['original_metadata']
        assert 'calibration' in dataset['original_metadata']

    def test_units_and_names(self, spectral_map_file):
        """Test proper units and axis names."""
        results = file_reader(spectral_map_file)
        dataset = results[0]
        
        # Check spatial axes
        assert dataset['axes'][0]['units'] == 'µm'
        assert dataset['axes'][1]['units'] == 'µm'
        
        # Check spectral axis
        assert dataset['axes'][2]['units'] == 'cm⁻¹'
        
        # Check navigation flags
        assert dataset['axes'][0]['navigate'] is True  # y
        assert dataset['axes'][1]['navigate'] is True  # x
        assert dataset['axes'][2]['navigate'] is False # spectral