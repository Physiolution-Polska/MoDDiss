#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import json

def Time_Intervals_read():
    """
    Read the intervals from json file. The intervals has to be ordered from low to high, include first and last point time
    """   

    input_file = open('time_interval.dat')
    input_str = input_file.read()
    Param = json.loads(input_str)
    list_of_periods = Param.values()
    list_of_periods.sort()
    
    return list_of_periods


