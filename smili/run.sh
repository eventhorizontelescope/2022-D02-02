#!/usr/bin/env bash
#
# Copyright (C) 2022 The Event Horizon Telescope Collaboration
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <https://www.gnu.org/licenses/>.


# Input uvfits directory of the data product release (e.g., 2022-XX-XXXX)
uvfitsdir="./uvfits"

# parameter table
parameter="./param_table/topset_deblur_2d_hops_3599_lo+hi_1parameter.csv"

#python $driver_script -p $ptabname -i ${inputdir}
./example_driver.py --uvfitsdir $uvfitsdir --parameter $parameter
