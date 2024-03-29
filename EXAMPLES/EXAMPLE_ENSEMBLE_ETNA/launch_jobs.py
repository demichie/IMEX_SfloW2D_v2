#!/usr/bin/env python
# -*- coding: utf-8 -*-

import os
import subprocess

dirs = next(os.walk('.'))[1]

print('Start computation')

command = "./run.sh"
processes = set()
max_processes = 5

for x in dirs:

    if "ensemble." in x:

        os.chdir(x)
        processes.add(subprocess.Popen([command, x]))
        if len(processes) >= max_processes:
            os.wait()
            processes.difference_update(
                [p for p in processes if p.poll() is not None])
        os.chdir("..")
