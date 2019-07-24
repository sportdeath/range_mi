#!/usr/bin/env python3

import matplotlib.pyplot as plt
import numpy as np

def main():
    xy = np.loadtxt("noisy_pdf.txt")
    x, y = xy[:,0], xy[:,1]

    plt.plot(x, y)
    plt.show()

if __name__ == "__main__":
    main()
