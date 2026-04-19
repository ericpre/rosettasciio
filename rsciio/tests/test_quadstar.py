# -*- coding: utf-8 -*-
# Copyright 2007-2026 The HyperSpy developers
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
from pathlib import Path

import numpy as np
import pytest

from rsciio.quadstar._api import (
    _build_datetime,
    _decode_bytes,
    _find_first_data_position,
    _read_general_header,
    _read_trace_headers,
    _read_trace_info,
    file_reader,
)

TEST_DATA_PATH = Path(__file__).parent / "data" / "quadstar"


class TestHelpers:
    def test_decode_bytes_normal(self):
        assert _decode_bytes(b"hello\x00\x00") == "hello"

    def test_decode_bytes_empty(self):
        assert _decode_bytes(b"\x00\x00") == ""

    def test_decode_bytes_non_bytes(self):
        assert _decode_bytes(42) == "42"

    def test_build_datetime_valid(self):
        header = {
            "year": 124,  # 2024
            "month": 3,
            "day": 15,
            "hour": 10,
            "minute": 30,
            "second": 45,
        }
        dt = _build_datetime(header)
        assert dt.year == 2024
        assert dt.month == 3
        assert dt.day == 15
        assert dt.hour == 10
        assert dt.minute == 30
        assert dt.second == 45

    def test_build_datetime_invalid(self):
        header = {
            "year": 0,
            "month": 0,
            "day": 0,
            "hour": 0,
            "minute": 0,
            "second": 0,
        }
        dt = _build_datetime(header)
        assert dt is None

    def test_find_first_data_position(self):
        headers = [
            {"type": 0x00, "data_position": 100},
            {"type": 0x11, "data_position": 500},
            {"type": 0x11, "data_position": 900},
        ]
        assert _find_first_data_position(headers) == 500

    def test_find_first_data_position_none(self):
        headers = [
            {"type": 0x00, "data_position": 100},
        ]
        assert _find_first_data_position(headers) is None


class TestReadTestSac:
    @pytest.fixture(autouse=True)
    def setup(self):
        self.filename = TEST_DATA_PATH / "test.sac"
        self.signals = file_reader(self.filename)

    def test_returns_list(self):
        assert isinstance(self.signals, list)
        assert len(self.signals) >= 1

    def test_signal_dict_keys(self):
        for sig in self.signals:
            assert "data" in sig
            assert "axes" in sig
            assert "metadata" in sig
            assert "original_metadata" in sig

    def test_data_is_numpy_array(self):
        for sig in self.signals:
            assert isinstance(sig["data"], np.ndarray)

    def test_data_is_float(self):
        for sig in self.signals:
            assert sig["data"].dtype == np.float32

    def test_axes_structure(self):
        for sig in self.signals:
            axes = sig["axes"]
            assert isinstance(axes, list)
            assert len(axes) >= 1
            for ax in axes:
                assert "name" in ax
                assert "units" in ax
                assert "size" in ax
                assert "index_in_array" in ax
                assert "navigate" in ax

    def test_signal_axis_has_calibration(self):
        for sig in self.signals:
            # The last axis should be the signal (M/Z) axis
            signal_ax = [a for a in sig["axes"] if not a["navigate"]]
            assert len(signal_ax) == 1
            ax = signal_ax[0]
            assert "offset" in ax
            assert "scale" in ax
            assert ax["scale"] > 0

    def test_metadata_general(self):
        for sig in self.signals:
            gen = sig["metadata"]["General"]
            assert "original_filename" in gen

    def test_metadata_signal(self, caplog):
        for sig in self.signals:
            with caplog.at_level(logging.WARNING):
                # Ignore not understood warnings about "signal_type" value
                assert sig["metadata"]["Signal"]["signal_type"] == "MS"

    def test_original_metadata_structure(self):
        for sig in self.signals:
            om = sig["original_metadata"]
            assert "general_header" in om
            assert "trace_info" in om
            assert "timestamps" in om

    def test_original_metadata_header_fields(self):
        om = self.signals[0]["original_metadata"]
        gh = om["general_header"]
        assert "n_timesteps" in gh
        assert "n_traces" in gh
        assert "timestep_length" in gh
        assert "username" in gh

    def test_original_metadata_trace_info_fields(self):
        om = self.signals[0]["original_metadata"]
        ti = om["trace_info"]
        assert "first_mass" in ti
        assert "scan_width" in ti
        assert "values_per_mass" in ti

    def test_data_shape_matches_axes(self):
        for sig in self.signals:
            data = sig["data"]
            axes = sig["axes"]
            if data.ndim == 1:
                assert len(axes) == 1
                assert axes[0]["size"] == data.shape[0]
            elif data.ndim == 2:
                assert len(axes) == 2
                for ax in axes:
                    assert ax["size"] == data.shape[ax["index_in_array"]]

    def test_no_all_nan_data(self):
        for sig in self.signals:
            assert not np.all(np.isnan(sig["data"]))


class TestReadAirdemoSac:
    @pytest.fixture(autouse=True)
    def setup(self):
        self.filename = TEST_DATA_PATH / "airdemo.sac"
        self.signals = file_reader(self.filename)

    def test_returns_list(self):
        assert isinstance(self.signals, list)
        assert len(self.signals) >= 1

    def test_signal_dict_keys(self):
        for sig in self.signals:
            assert "data" in sig
            assert "axes" in sig
            assert "metadata" in sig
            assert "original_metadata" in sig

    def test_data_is_numpy_array(self):
        for sig in self.signals:
            assert isinstance(sig["data"], np.ndarray)
            assert sig["data"].dtype == np.float32

    def test_data_shape_matches_axes(self):
        for sig in self.signals:
            data = sig["data"]
            for ax in sig["axes"]:
                assert ax["size"] == data.shape[ax["index_in_array"]]

    def test_mass_axis_positive(self):
        for sig in self.signals:
            signal_ax = [a for a in sig["axes"] if not a["navigate"]]
            for ax in signal_ax:
                assert ax["offset"] >= 0
                assert ax["scale"] > 0

    def test_timestamps_in_original_metadata(self):
        for sig in self.signals:
            ts = sig["original_metadata"]["timestamps"]
            assert isinstance(ts, list)
            assert len(ts) >= 1
            # Timestamps should be monotonically non-decreasing
            for i in range(1, len(ts)):
                assert ts[i] >= ts[i - 1]


class TestReadGeneralHeader:
    def test_general_header(self):
        with open(TEST_DATA_PATH / "test.sac", "rb") as f:
            buf = f.read()
        gen = _read_general_header(buf)
        assert isinstance(gen, dict)
        assert gen["n_timesteps"] > 0
        assert gen["n_traces"] > 0
        assert gen["timestep_length"] > 0

    def test_trace_headers(self):
        with open(TEST_DATA_PATH / "test.sac", "rb") as f:
            buf = f.read()
        gen = _read_general_header(buf)
        headers = _read_trace_headers(buf, gen["n_traces"])
        assert len(headers) == gen["n_traces"]
        # At least one trace should be type 0x11 (ScanAnalog)
        types = [h["type"] for h in headers]
        assert 0x11 in types

    def test_trace_info(self):
        with open(TEST_DATA_PATH / "test.sac", "rb") as f:
            buf = f.read()
        gen = _read_general_header(buf)
        headers = _read_trace_headers(buf, gen["n_traces"])
        for h in headers:
            if h["type"] == 0x11:
                info = _read_trace_info(buf, h["info_position"])
                assert info["scan_width"] > 0
                assert info["values_per_mass"] > 0
                assert info["first_mass"] >= 0
                break


class TestLazyNotSupported:
    def test_lazy_raises(self):
        with pytest.raises(NotImplementedError, match="Lazy loading"):
            file_reader(TEST_DATA_PATH / "test.sac", lazy=True)


class TestInvalidFile:
    def test_empty_file(self, tmp_path):
        empty = tmp_path / "empty.sac"
        empty.write_bytes(b"")
        with pytest.raises(Exception):
            file_reader(empty)

    def test_too_small_file(self, tmp_path):
        small = tmp_path / "small.sac"
        small.write_bytes(b"\x00" * 10)
        with pytest.raises(Exception):
            file_reader(small)
