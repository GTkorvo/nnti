#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
#  Copyright (c) 2014, Nico Schlömer, <nico.schloemer@gmail.com>
#  All rights reserved.
#
#  This file is part of Sundance.
#
#  Redistribution and use in source and binary forms, with or without
#  modification, are permitted provided that the following conditions are met:
#
#  1. Redistributions of source code must retain the above copyright notice,
#  this list of conditions and the following disclaimer.
#
#  2. Redistributions in binary form must reproduce the above copyright notice,
#  this list of conditions and the following disclaimer in the documentation
#  and/or other materials provided with the distribution.
#
#  3. Neither the name of the copyright holder nor the names of its
#  contributors may be used to endorse or promote products derived from this
#  software without specific prior written permission.
#
#  THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
#  AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
#  IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
#  ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE
#  LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
#  CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
#  SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
#  INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
#  CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
#  ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
#  POSSIBILITY OF SUCH DAMAGE.
#

from mpl_toolkits.mplot3d import Axes3D
from matplotlib import cm
#from matplotlib.ticker import LinearLocator, FormatStrFormatter
import matplotlib.pyplot as plt
import numpy as np


def rosen(X, Y, alpha=1.0):
    """The Rosenbrock function"""
    return alpha * (Y - X**2)**2 + (X - 1.0)**2


def plot_rosen(alpha):
    fig = plt.figure()
    ax = fig.gca(projection='3d')
    ax.set_aspect('equal')
    x = np.arange(-2, 2, 5.0e-2)
    y = np.arange(0, 3, 5.0e-2)
    X, Y = np.meshgrid(x, y)
    Z = rosen(X, Y, alpha=alpha)
    ax.plot_surface(
        X, Y, Z,
        rstride=1,
        cstride=1,
        cmap=cm.jet_r,
        linewidth=0,
        antialiased=False
        )


if __name__ == '__main__':
    plot_rosen(alpha=1.0)
    plot_rosen(alpha=10.0)
    plot_rosen(alpha=100.0)
    plt.show()
