#!/usr/bin/env python3.8
import numpy as np
import matplotlib.pyplot as plt

def plot_sa_bars(ax, main_effect, total_effect, names):  # ax required for dycoms multiplot
    plt.subplots_adjust(left=0.15,right=0.8)
    
    interactions = total_effect - main_effect

    bars = {
        "Main effect": [main_effect*100, "black"],
        "Interactions": [interactions*100, "white"]
    }

    
    bottom=np.zeros(len(main_effect))
    for label, val in bars.items():
        p = ax.bar(names, val[0], label=label, bottom=bottom, color=val[1], edgecolor = "black")
        bottom += val[0]

    ax.set_ylabel("Percentage \nof variance (%)")   #turned off for dycoms multiplot
    ax.set_ylim((0,100))
    ax.legend(loc=(1.01, 0), handlelength=0.75, frameon=False, fontsize=15,handletextpad=0.5)   #turned off for dycoms multiplot

    for n, v in zip(names, main_effect):
        print(f"{n} has average main effect contribution of {v*100:0.2f}%")
    return ax

def plot_validation(ax, loo_prediction, actualval, max_line, title, font_size):
    mean = np.array(loo_prediction[0])
    sd = np.array(loo_prediction[1])
    lower95 = mean - 1.96*sd
    upper95 = mean + 1.96*sd

    minX = min(actualval)
    maxX = max(actualval)
    minY = min(lower95)
    maxY = max(upper95)
    minXY = min(minY, minX)
    maxXY = max(maxY, maxX)

    seq = np.linspace(minXY,max_line,10)
    errors = []
    colours = []
    count_fail = 0
    for v in range(len(upper95)):
        errors.append((upper95[v] - lower95[v])/2)
        if actualval[v] > lower95[v] and actualval[v] < upper95[v]:
            colours.append("black")
        else:
            colours.append("red")
            count_fail+=1

    per_diff = np.mean([e*100/m for e,m in zip(errors,mean)])
    # print(per_diff)

    line_of_equality, = ax.plot(seq, seq, c='black',label='line of equality',linestyle='--')
    points = ax.scatter(actualval, mean, color=colours,s=50,marker='o')

    for val, error, color, l95, u95 in zip(actualval, errors, colours, lower95, upper95):
        line, = ax.plot((val, val),(l95, u95),c=color, linewidth=1)

    rmse = np.mean((actualval - mean)**2)**0.5
    r_squared = calc_r_squared(actualval, mean)
    ax.text(0.03,0.9,'NRMSE = {:0.2f}'.format(rmse/np.mean(mean)), transform=ax.transAxes, fontsize=font_size, horizontalalignment="left", verticalalignment="bottom")
    ax.text(0.03,0.81,'Pass = {:2.0f}%'.format(100 - (count_fail*100/len(actualval))), transform=ax.transAxes, fontsize=font_size, horizontalalignment="left", verticalalignment="bottom")
    ax.text(0.03,0.72,'r$^2$ = {:0.2f}'.format(r_squared), transform=ax.transAxes, fontsize=font_size, horizontalalignment="left", verticalalignment="bottom")
    ax.set_title(title)
    return rmse

def calc_r_squared(model_values, emulator_values):
    em_mean = np.mean(emulator_values)
    residuals = []
    squares = []
    for m_val, e_val in zip(model_values, emulator_values):
        residuals.append((e_val - m_val)**2)
        squares.append((e_val - em_mean)**2)
    r_squared = 1 - (np.sum(residuals))/(np.sum(squares))
    return r_squared

    print(len(mean))
    print(len(actualval))
    rmse = np.mean((actualval - mean)**2)**0.5
    ax.text(0.03,0.9,'Pass = {:2.0f}%'.format(100 - (count_fail*100/len(actualval))), transform=ax.transAxes, fontsize=font_size, horizontalalignment="left", verticalalignment="bottom")
    ax.set_title(title)
    return rmse