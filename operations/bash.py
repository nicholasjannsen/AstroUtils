#!/usr/bin/env python3

import sys
import time


def compilation(i, i_max, text):
    """
    This function print out a compilation time menu bar in the terminal
    """
    percent = (i + 1) / (i_max * 1.0) * 100
    # print int(percent/2), 50-int(percent/2)
    # We here divide by 2 as the length of the bar is only 50 characters:
    bar = "[" + "-" * int(percent/2) + '>' + " " *(50 - int(percent/2)) + "] {}% {}".format(int(percent), text)
    sys.stdout.write(u"\u001b[1000D" +  bar)
    sys.stdout.flush()


def loading(i, i_max):
    """
    This function print the loading status
    """
    sys.stdout.write(u"\u001b[1000D")
    sys.stdout.flush()
    time.sleep(1)
    sys.stdout.write(str(i + 1) + "Loding... %")
    sys.stdout.flush()
