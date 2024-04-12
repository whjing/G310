# !/user/bin/env python
# -*- coding: utf-8 -*-
import sys
import os
res = "/Users/jing/academic/module"
#res = os.path.dirname(os.path.dirname(__file__))
sys.path.append(res)
print(res)
print('Module path successfully appended.')