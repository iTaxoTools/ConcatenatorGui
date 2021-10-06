#!/usr/bin/env python3

from random import choices, randint

import sys


def fetch(list, index, default):
    try:
        return list[index]
    except IndexError:
        return default


def sequence(low=10, high=10):
    return ''.join(choices('AGCT-', k=randint(low, high)))


def main():
    sample_count = int(fetch(sys.argv, 1, randint(10, 20)))
    charset_count = int(fetch(sys.argv, 2, randint(1, 10)))
    sample_prefix = fetch(sys.argv, 3, 'sample')
    charset_prefix = fetch(sys.argv, 4, 'charset')

    print('species\tspecimen-voucher\tlocality', end='')
    for charset in range(charset_count):
        print(f'\tsequence_{charset_prefix}_{str(charset)}', end='')
    print()
    for sample in range(sample_count):
        print(f'{sample_prefix}_{str(sample)}', end='')
        print(f'\tspecimen_{str(sample)}', end='')
        print(f'\tlocality_{str(sample)}', end='')
        for charset in range(charset_count):
            print(f'\t{sequence()}', end='')
        print()


if __name__ == '__main__':
    main()
