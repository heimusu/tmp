# -*- coding:utf-8 -*-
import cmath

if __name__ == "__main__":
    a = abs(1 + 2j)**2
    b = 11
    c = 0
    d = cmath.pi
    print a
    for i in range(0,b-1,2):
	print i
	c += i
    print c
    print d
