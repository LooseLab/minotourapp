/**
 * minoTour Highcharts theme — Obsidian & Ether (see design.md).
 * Replaces stock Grid-light: Manrope matches template_base; no extra webfont load.
 */
'use strict';

/* global Highcharts */

Highcharts.theme = {
  colors: [
    '#004B44',
    '#0d9488',
    '#14b8a6',
    '#059669',
    '#0891b2',
    '#475569',
    '#ca8a04',
    '#0ea5e9',
    '#8b5cf6'
  ],
  chart: {
    backgroundColor: 'transparent',
    borderWidth: 0,
    plotBorderWidth: 0,
    style: {
      fontFamily: '"Manrope", -apple-system, BlinkMacSystemFont, "Segoe UI", sans-serif',
      fontSize: '12px',
      color: '#334155'
    }
  },
  title: {
    align: 'left',
    margin: 12,
    style: {
      color: '#00322d',
      fontSize: '15px',
      fontWeight: '600',
      fontFamily: '"Manrope", -apple-system, BlinkMacSystemFont, "Segoe UI", sans-serif',
      textTransform: 'none'
    }
  },
  subtitle: {
    style: {
      color: '#64748b',
      fontSize: '12px',
      fontWeight: '400'
    }
  },
  xAxis: {
    lineColor: '#e2e8f0',
    tickColor: '#e2e8f0',
    gridLineColor: '#f1f5f9',
    gridLineWidth: 1,
    labels: {
      style: {
        color: '#64748b',
        fontSize: '11px',
        fontWeight: '500'
      }
    },
    title: {
      style: {
        color: '#475569',
        fontSize: '11px',
        fontWeight: '600',
        textTransform: 'none'
      }
    }
  },
  yAxis: {
    lineColor: '#e2e8f0',
    tickColor: '#e2e8f0',
    gridLineColor: '#f1f5f9',
    gridLineWidth: 1,
    labels: {
      style: {
        color: '#64748b',
        fontSize: '11px',
        fontWeight: '500'
      }
    },
    title: {
      style: {
        color: '#475569',
        fontSize: '11px',
        fontWeight: '600',
        textTransform: 'none'
      }
    }
  },
  legend: {
    itemStyle: {
      color: '#475569',
      fontWeight: '500',
      fontSize: '12px'
    },
    itemHoverStyle: {
      color: '#00322d'
    },
    itemHiddenStyle: {
      color: '#94a3b8',
      textDecoration: 'line-through'
    }
  },
  tooltip: {
    backgroundColor: '#ffffff',
    borderColor: '#e2e8f0',
    borderRadius: 6,
    borderWidth: 1,
    shadow: false,
    style: {
      color: '#1e293b',
      fontSize: '12px'
    }
  },
  plotOptions: {
    series: {
      borderWidth: 0,
      lineWidth: 2,
      states: {
        hover: {
          lineWidth: 2,
          halo: {
            size: 6,
            opacity: 0.12
          }
        }
      }
    },
    column: {
      borderWidth: 0,
      borderRadius: 3,
      groupPadding: 0.08
    },
    bar: {
      borderWidth: 0,
      borderRadius: 3
    },
    pie: {
      borderWidth: 0,
      dataLabels: {
        style: {
          fontWeight: '500',
          color: '#475569'
        }
      }
    },
    boxplot: {
      lineWidth: 1,
      medianWidth: 2,
      stemWidth: 1,
      whiskerWidth: 1
    },
    candlestick: {
      lineColor: '#475569'
    }
  },
  credits: {
    enabled: false
  },
  navigation: {
    buttonOptions: {
      theme: {
        fill: '#f8fafc',
        stroke: '#e2e8f0',
        style: {
          color: '#475569'
        },
        states: {
          hover: {
            fill: '#f1f5f9',
            stroke: '#cbd5e1'
          },
          select: {
            fill: '#ecfdf5',
            stroke: '#004B44',
            style: {
              color: '#00322d'
            }
          }
        }
      }
    }
  },
  rangeSelector: {
    buttonTheme: {
      fill: '#f8fafc',
      stroke: '#e2e8f0',
      style: {
        color: '#475569',
        fontWeight: '500'
      },
      states: {
        hover: {
          fill: '#ecfdf5',
          stroke: '#99f6e4'
        },
        select: {
          fill: '#004B44',
          stroke: '#004B44',
          style: {
            color: '#ffffff',
            fontWeight: '600'
          }
        }
      },
      width: 40,
      height: 22,
      padding: 4,
      r: 4
    },
    inputBoxBorderColor: '#e2e8f0',
    inputBoxBackgroundColor: '#ffffff',
    inputStyle: {
      color: '#334155',
      fontWeight: '500'
    },
    labelStyle: {
      color: '#64748b'
    }
  },
  navigator: {
    outlineColor: '#e2e8f0',
    maskFill: 'rgba(0, 75, 68, 0.08)',
    series: {
      color: '#94a3b8',
      lineColor: '#64748b'
    },
    xAxis: {
      gridLineColor: '#f1f5f9'
    }
  },
  scrollbar: {
    barBackgroundColor: '#f1f5f9',
    barBorderColor: '#e2e8f0',
    buttonBackgroundColor: '#f8fafc',
    buttonBorderColor: '#e2e8f0',
    rifleColor: '#94a3b8',
    trackBackgroundColor: '#f8fafc',
    trackBorderColor: '#e2e8f0'
  },
  background2: '#f1f5f9',
  contrastTextColor: '#00322d'
};

Highcharts.setOptions(Highcharts.theme);
/* Snapshot for toggling back from dark mode (minotour-theme.js). */
if (typeof Highcharts !== 'undefined' && Highcharts.theme) {
  try {
    window.__minotourHcLightTheme = Highcharts.merge(true, {}, Highcharts.theme);
  } catch (e) {
    window.__minotourHcLightTheme = Highcharts.theme;
  }
}
