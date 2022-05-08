#!/usr/bin/env python3

from colorama import Fore, Style



def warning(message):
    print(Style.BRIGHT + Fore.YELLOW + message + Style.RESET_ALL)


def error(message):
    print(Style.BRIGHT + Fore.RED + message + Style.RESET_ALL)
    exit()
