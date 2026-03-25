/**
 * Site appearance: html[data-theme="light|dark"], localStorage minotour-theme.
 * Depends on Highcharts + grid-light-custom for __minotourHcLightTheme snapshot.
 */
(function (window) {
  'use strict'

  var STORAGE_KEY = 'minotour-theme'
  var VALID_THEMES = ['light', 'dark', 'hacienda', 'sf1950s', 'paris1940s', 'hailmary', 'stoneage', 'disney', 'pixar', 'telextext', 'flowerpower1960s', 'marblegold', 'rusticwood', 'scifi1950s', 'drwho', 'playschool', 'tron', 'starwarsanh', 'academic', 'mutedscientist', 'terminal', 'kimiko', 'londoncalling', 'grumpychemist', 'oxford', 'cambridge']

  /** Themes with dark chrome; enables Tailwind dark: variants (html.dark) for settings drawer etc. Keep in sync with template_base.html inline script. */
  var DARK_UI_THEMES = {
    dark: true,
    hacienda: true,
    hailmary: true,
    telextext: true,
    scifi1950s: true,
    drwho: true,
    tron: true,
    starwarsanh: true,
    terminal: true,
    grumpychemist: true,
    londoncalling: true
  }

  var HIGHCHARTS_DARK = {
    colors: [
      '#5eead4',
      '#34d399',
      '#2dd4bf',
      '#6ee7b7',
      '#38bdf8',
      '#a78bfa',
      '#fbbf24',
      '#fb7185',
      '#94a3b8'
    ],
    chart: {
      backgroundColor: 'transparent',
      borderWidth: 0,
      plotBorderWidth: 0,
      style: {
        fontFamily: '"IBM Plex Sans", -apple-system, BlinkMacSystemFont, "Segoe UI", sans-serif',
        fontSize: '12px',
        color: '#cbd5e1'
      }
    },
    title: {
      align: 'left',
      margin: 12,
      style: {
        color: '#ecfdf5',
        fontSize: '15px',
        fontWeight: '600',
        fontFamily: '"IBM Plex Sans", -apple-system, BlinkMacSystemFont, "Segoe UI", sans-serif',
        textTransform: 'none'
      }
    },
    subtitle: {
      style: {
        color: '#94a3b8',
        fontSize: '12px',
        fontWeight: '400'
      }
    },
    xAxis: {
      lineColor: '#334155',
      tickColor: '#334155',
      gridLineColor: '#1e293b',
      gridLineWidth: 1,
      labels: {
        style: {
          color: '#94a3b8',
          fontSize: '11px',
          fontWeight: '500'
        }
      },
      title: {
        style: {
          color: '#cbd5e1',
          fontSize: '11px',
          fontWeight: '600',
          textTransform: 'none'
        }
      }
    },
    yAxis: {
      lineColor: '#334155',
      tickColor: '#334155',
      gridLineColor: '#1e293b',
      gridLineWidth: 1,
      labels: {
        style: {
          color: '#94a3b8',
          fontSize: '11px',
          fontWeight: '500'
        }
      },
      title: {
        style: {
          color: '#cbd5e1',
          fontSize: '11px',
          fontWeight: '600',
          textTransform: 'none'
        }
      }
    },
    legend: {
      itemStyle: {
        color: '#cbd5e1',
        fontWeight: '500',
        fontSize: '12px'
      },
      itemHoverStyle: {
        color: '#ecfdf5'
      },
      itemHiddenStyle: {
        color: '#64748b',
        textDecoration: 'line-through'
      }
    },
    tooltip: {
      backgroundColor: '#1e293b',
      borderColor: '#334155',
      borderRadius: 6,
      borderWidth: 1,
      shadow: false,
      style: {
        color: '#f1f5f9',
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
              opacity: 0.15
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
            color: '#cbd5e1'
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
        lineColor: '#94a3b8'
      }
    },
    credits: {
      enabled: false
    },
    navigation: {
      buttonOptions: {
        theme: {
          fill: '#1e293b',
          stroke: '#334155',
          style: {
            color: '#e2e8f0'
          },
          states: {
            hover: {
              fill: '#334155',
              stroke: '#475569'
            },
            select: {
              fill: '#0f766e',
              stroke: '#14b8a6',
              style: {
                color: '#ecfdf5'
              }
            }
          }
        }
      }
    },
    rangeSelector: {
      buttonTheme: {
        fill: '#1e293b',
        stroke: '#334155',
        style: {
          color: '#cbd5e1',
          fontWeight: '500'
        },
        states: {
          hover: {
            fill: '#334155',
            stroke: '#475569'
          },
          select: {
            fill: '#0f766e',
            stroke: '#14b8a6',
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
      inputBoxBorderColor: '#334155',
      inputBoxBackgroundColor: '#0f172a',
      inputStyle: {
        color: '#e2e8f0',
        fontWeight: '500'
      },
      labelStyle: {
        color: '#94a3b8'
      }
    },
    navigator: {
      outlineColor: '#334155',
      maskFill: 'rgba(15, 23, 42, 0.45)',
      series: {
        color: '#64748b',
        lineColor: '#475569'
      },
      xAxis: {
        gridLineColor: '#1e293b'
      }
    },
    scrollbar: {
      barBackgroundColor: '#1e293b',
      barBorderColor: '#334155',
      buttonBackgroundColor: '#0f172a',
      buttonBorderColor: '#334155',
      rifleColor: '#64748b',
      trackBackgroundColor: '#0f172a',
      trackBorderColor: '#334155'
    },
    background2: '#334155',
    contrastTextColor: '#ecfdf5'
  }

  var HIGHCHARTS_HACIENDA = {
    colors: [
      '#ffd60a',
      '#ff006e',
      '#44d7ff',
      '#f3ff4a',
      '#ff7f11',
      '#9d4edd',
      '#39ffb6',
      '#ff5c8a',
      '#ffe066'
    ],
    chart: {
      backgroundColor: 'transparent',
      borderWidth: 0,
      plotBorderWidth: 0,
      style: {
        fontFamily: '"Barlow Condensed", "IBM Plex Sans", -apple-system, BlinkMacSystemFont, "Segoe UI", sans-serif',
        fontSize: '12px',
        color: '#ffeef8'
      }
    },
    title: {
      align: 'left',
      margin: 12,
      style: {
        color: '#ffd60a',
        fontSize: '15px',
        fontWeight: '700'
      }
    },
    subtitle: {
      style: {
        color: '#ffb3d7',
        fontSize: '12px'
      }
    },
    xAxis: {
      lineColor: '#ffd60a',
      tickColor: '#ffd60a',
      gridLineColor: 'rgba(255, 214, 10, 0.22)',
      labels: { style: { color: '#ffcae5' } },
      title: { style: { color: '#ffe6f3' } }
    },
    yAxis: {
      lineColor: '#ffd60a',
      tickColor: '#ffd60a',
      gridLineColor: 'rgba(255, 214, 10, 0.22)',
      labels: { style: { color: '#ffcae5' } },
      title: { style: { color: '#ffe6f3' } }
    },
    legend: {
      itemStyle: { color: '#ffd7ee', fontWeight: '600' },
      itemHoverStyle: { color: '#fff46b' },
      itemHiddenStyle: { color: '#8c4f90' }
    },
    tooltip: {
      backgroundColor: '#140816',
      borderColor: '#ffd60a',
      borderRadius: 6,
      borderWidth: 1,
      shadow: false,
      style: { color: '#fff0f8', fontSize: '12px' }
    },
    plotOptions: {
      series: { borderWidth: 0, lineWidth: 2 },
      column: { borderWidth: 0, borderRadius: 3, groupPadding: 0.08 },
      bar: { borderWidth: 0, borderRadius: 3 },
      pie: { borderWidth: 0, dataLabels: { style: { color: '#ffd7ee' } } }
    },
    credits: { enabled: false }
  }

  var HIGHCHARTS_SF1950S = {
    colors: [
      '#d45c3f',
      '#4f7a78',
      '#f2c572',
      '#9d5e4f',
      '#8e9a6e',
      '#b66a5e',
      '#5b6f8e',
      '#c7874c',
      '#7d8a5c'
    ],
    chart: {
      backgroundColor: 'transparent',
      borderWidth: 0,
      plotBorderWidth: 0,
      style: {
        fontFamily: '"Source Sans 3", "Manrope", -apple-system, BlinkMacSystemFont, "Segoe UI", sans-serif',
        fontSize: '12px',
        color: '#2d2a26'
      }
    },
    title: {
      align: 'left',
      margin: 12,
      style: {
        color: '#5a3527',
        fontSize: '15px',
        fontWeight: '700'
      }
    },
    subtitle: {
      style: {
        color: '#6f6256',
        fontSize: '12px'
      }
    },
    xAxis: {
      lineColor: '#b9a78e',
      tickColor: '#b9a78e',
      gridLineColor: 'rgba(114, 132, 101, 0.22)',
      labels: { style: { color: '#5b4d3f' } },
      title: { style: { color: '#5a3527' } }
    },
    yAxis: {
      lineColor: '#b9a78e',
      tickColor: '#b9a78e',
      gridLineColor: 'rgba(114, 132, 101, 0.22)',
      labels: { style: { color: '#5b4d3f' } },
      title: { style: { color: '#5a3527' } }
    },
    legend: {
      itemStyle: { color: '#4a443d', fontWeight: '600' },
      itemHoverStyle: { color: '#5a3527' },
      itemHiddenStyle: { color: '#8f857a' }
    },
    tooltip: {
      backgroundColor: '#f6efe2',
      borderColor: '#c8b79f',
      borderRadius: 6,
      borderWidth: 1,
      shadow: false,
      style: { color: '#2e2a25', fontSize: '12px' }
    },
    plotOptions: {
      series: { borderWidth: 0, lineWidth: 2 },
      column: { borderWidth: 0, borderRadius: 3, groupPadding: 0.08 },
      bar: { borderWidth: 0, borderRadius: 3 }
    },
    credits: { enabled: false }
  }

  var HIGHCHARTS_PARIS1940S = {
    colors: ['#7a4f2f', '#4e5b53', '#9f7b4a', '#5f4733', '#6e6b5a', '#7b6b4d', '#8c5a3c'],
    chart: {
      backgroundColor: 'transparent',
      style: {
        color: '#2f2923',
        fontFamily: '"Libre Baskerville", "Source Sans 3", Georgia, serif'
      }
    },
    title: {
      style: {
        color: '#3f2f22',
        fontWeight: '700',
        fontFamily: '"EB Garamond", "Cormorant Garamond", Georgia, serif'
      }
    },
    subtitle: { style: { color: '#6a5b4c', fontFamily: '"Libre Baskerville", "Source Sans 3", Georgia, serif' } },
    xAxis: { lineColor: '#b3a48f', tickColor: '#b3a48f', gridLineColor: 'rgba(122, 79, 47, 0.16)', labels: { style: { color: '#5c5043' } } },
    yAxis: { lineColor: '#b3a48f', tickColor: '#b3a48f', gridLineColor: 'rgba(122, 79, 47, 0.16)', labels: { style: { color: '#5c5043' } } },
    legend: { itemStyle: { color: '#4a4036' } },
    tooltip: { backgroundColor: '#f2e6d2', borderColor: '#c8b79f', style: { color: '#2f2923' } },
    credits: { enabled: false }
  }

  var HIGHCHARTS_HAILMARY = {
    colors: ['#7df9ff', '#20d8ff', '#56f000', '#f4ff61', '#68a8ff', '#b56dff', '#40f3d5'],
    chart: {
      backgroundColor: 'transparent',
      style: {
        color: '#d9f6ff',
        fontFamily: '"Exo 2", "IBM Plex Sans", -apple-system, BlinkMacSystemFont, "Segoe UI", sans-serif'
      }
    },
    title: {
      style: {
        color: '#c2ff7a',
        fontWeight: '700',
        fontFamily: '"Orbitron", "Exo 2", "IBM Plex Sans", -apple-system, BlinkMacSystemFont, "Segoe UI", sans-serif'
      }
    },
    subtitle: { style: { color: '#88d6ff', fontFamily: '"Exo 2", "IBM Plex Sans", -apple-system, BlinkMacSystemFont, "Segoe UI", sans-serif' } },
    xAxis: { lineColor: '#2b4f6a', tickColor: '#2b4f6a', gridLineColor: 'rgba(125, 249, 255, 0.2)', labels: { style: { color: '#8fdfff' } } },
    yAxis: { lineColor: '#2b4f6a', tickColor: '#2b4f6a', gridLineColor: 'rgba(125, 249, 255, 0.2)', labels: { style: { color: '#8fdfff' } } },
    legend: { itemStyle: { color: '#b8ecff' } },
    tooltip: { backgroundColor: '#0b1322', borderColor: '#20d8ff', style: { color: '#e7fbff' } },
    credits: { enabled: false }
  }

  var HIGHCHARTS_STONEAGE = {
    colors: ['#b3532e', '#7f5f3f', '#d18a3f', '#5f7f4f', '#9c6b36', '#c13f2f', '#7a4f2f'],
    chart: {
      backgroundColor: 'transparent',
      style: {
        color: '#2b2218',
        fontFamily: '"Source Sans 3", "Manrope", -apple-system, BlinkMacSystemFont, "Segoe UI", sans-serif'
      }
    },
    title: {
      style: {
        color: '#4b2f1c',
        fontWeight: '700',
        fontFamily: '"Bangers", "Archivo Black", "Barlow Condensed", sans-serif'
      }
    },
    subtitle: { style: { color: '#6d5845' } },
    xAxis: { lineColor: '#b99873', tickColor: '#b99873', gridLineColor: 'rgba(75, 47, 28, 0.18)', labels: { style: { color: '#5c4a39' } } },
    yAxis: { lineColor: '#b99873', tickColor: '#b99873', gridLineColor: 'rgba(75, 47, 28, 0.18)', labels: { style: { color: '#5c4a39' } } },
    legend: { itemStyle: { color: '#4b3a2a' } },
    tooltip: { backgroundColor: '#f5e4c4', borderColor: '#c79c64', style: { color: '#2b2218' } },
    credits: { enabled: false }
  }

  var HIGHCHARTS_DISNEY = {
    colors: ['#1e6fd9', '#56b0ff', '#ffd166', '#ff7aa2', '#70d6ff', '#7bc96f', '#8a7dff'],
    chart: {
      backgroundColor: 'transparent',
      style: {
        color: '#173055',
        fontFamily: '"Nunito", "Manrope", -apple-system, BlinkMacSystemFont, "Segoe UI", sans-serif'
      }
    },
    title: {
      style: {
        color: '#1e4f9e',
        fontWeight: '700',
        fontFamily: '"Baloo 2", "Nunito", "Manrope", sans-serif'
      }
    },
    subtitle: { style: { color: '#3f6ea8' } },
    xAxis: { lineColor: '#9fc8f8', tickColor: '#9fc8f8', gridLineColor: 'rgba(30, 111, 217, 0.18)', labels: { style: { color: '#3f6ea8' } } },
    yAxis: { lineColor: '#9fc8f8', tickColor: '#9fc8f8', gridLineColor: 'rgba(30, 111, 217, 0.18)', labels: { style: { color: '#3f6ea8' } } },
    legend: { itemStyle: { color: '#2a4e86' } },
    tooltip: { backgroundColor: '#e8f4ff', borderColor: '#9fc8f8', style: { color: '#173055' } },
    credits: { enabled: false }
  }

  var HIGHCHARTS_PIXAR = {
    colors: ['#0b5fff', '#00b3ff', '#3ad4c4', '#ffb100', '#ff6f61', '#6f7dff', '#7cc36f'],
    chart: {
      backgroundColor: 'transparent',
      style: {
        color: '#10243f',
        fontFamily: '"Nunito", "Manrope", -apple-system, BlinkMacSystemFont, "Segoe UI", sans-serif'
      }
    },
    title: {
      style: {
        color: '#0b5fff',
        fontWeight: '800',
        fontFamily: '"Nunito", "Manrope", -apple-system, BlinkMacSystemFont, "Segoe UI", sans-serif'
      }
    },
    subtitle: { style: { color: '#3f6ea8' } },
    xAxis: { lineColor: '#8bc5ff', tickColor: '#8bc5ff', gridLineColor: 'rgba(11, 95, 255, 0.16)', labels: { style: { color: '#2f4f78' } } },
    yAxis: { lineColor: '#8bc5ff', tickColor: '#8bc5ff', gridLineColor: 'rgba(11, 95, 255, 0.16)', labels: { style: { color: '#2f4f78' } } },
    legend: { itemStyle: { color: '#1f3861' } },
    tooltip: { backgroundColor: '#eef7ff', borderColor: '#8bc5ff', style: { color: '#10243f' } },
    credits: { enabled: false }
  }

  var HIGHCHARTS_TELEXTEXT = {
    colors: ['#00ffff', '#ffff00', '#ff00ff', '#00ff00', '#ff3030', '#40a0ff', '#ffffff'],
    chart: {
      backgroundColor: 'transparent',
      style: {
        color: '#e8faff',
        fontFamily: '"VT323", "IBM Plex Sans", monospace'
      }
    },
    title: {
      style: {
        color: '#ffff00',
        fontWeight: '700',
        fontFamily: '"VT323", "IBM Plex Sans", monospace'
      }
    },
    subtitle: { style: { color: '#8adfff' } },
    xAxis: { lineColor: '#00ffff', tickColor: '#00ffff', gridLineColor: 'rgba(0, 255, 255, 0.2)', labels: { style: { color: '#a7eeff' } } },
    yAxis: { lineColor: '#00ffff', tickColor: '#00ffff', gridLineColor: 'rgba(0, 255, 255, 0.2)', labels: { style: { color: '#a7eeff' } } },
    legend: { itemStyle: { color: '#d5f7ff' } },
    tooltip: { backgroundColor: '#05101c', borderColor: '#00ffff', style: { color: '#e8faff' } },
    credits: { enabled: false }
  }

  var HIGHCHARTS_FLOWERPOWER1960S = {
    colors: ['#ff4fa1', '#ff9f1c', '#ffd60a', '#7ae582', '#47b8ff', '#9b5de5', '#ff6b6b'],
    chart: {
      backgroundColor: 'transparent',
      style: {
        color: '#3b235a',
        fontFamily: '"Nunito", "Manrope", -apple-system, BlinkMacSystemFont, "Segoe UI", sans-serif'
      }
    },
    title: {
      style: {
        color: '#9b2fc9',
        fontWeight: '700',
        fontFamily: '"Pacifico", "Baloo 2", "Nunito", sans-serif'
      }
    },
    subtitle: { style: { color: '#5d3a8f' } },
    xAxis: { lineColor: '#c58cff', tickColor: '#c58cff', gridLineColor: 'rgba(155, 47, 201, 0.16)', labels: { style: { color: '#5d3a8f' } } },
    yAxis: { lineColor: '#c58cff', tickColor: '#c58cff', gridLineColor: 'rgba(155, 47, 201, 0.16)', labels: { style: { color: '#5d3a8f' } } },
    legend: { itemStyle: { color: '#4a2f70' } },
    tooltip: { backgroundColor: '#fff2ff', borderColor: '#d77bff', style: { color: '#3b235a' } },
    credits: { enabled: false }
  }

  var HIGHCHARTS_MARBLEGOLD = {
    colors: ['#c8a24a', '#9f7c2c', '#d8c08a', '#6f6f73', '#bfa56a', '#8e8e94', '#d4af37'],
    chart: {
      backgroundColor: 'transparent',
      style: {
        color: '#2f3136',
        fontFamily: '"Lora", "Source Sans 3", Georgia, serif'
      }
    },
    title: {
      style: {
        color: '#8a6a23',
        fontWeight: '700',
        fontFamily: '"Cinzel", "Lora", Georgia, serif'
      }
    },
    subtitle: { style: { color: '#6b6560' } },
    xAxis: { lineColor: '#c8b27a', tickColor: '#c8b27a', gridLineColor: 'rgba(159, 124, 44, 0.16)', labels: { style: { color: '#6b6560' } } },
    yAxis: { lineColor: '#c8b27a', tickColor: '#c8b27a', gridLineColor: 'rgba(159, 124, 44, 0.16)', labels: { style: { color: '#6b6560' } } },
    legend: { itemStyle: { color: '#4d4f55' } },
    tooltip: { backgroundColor: '#f7f6f2', borderColor: '#c8b27a', style: { color: '#2f3136' } },
    credits: { enabled: false }
  }

  var HIGHCHARTS_RUSTICWOOD = {
    colors: ['#8c5a3c', '#a97142', '#6f8a52', '#c47b43', '#7a4f2f', '#9e6a3a', '#4f6a3e'],
    chart: {
      backgroundColor: 'transparent',
      style: {
        color: '#2e2218',
        fontFamily: '"Lora", "Source Sans 3", Georgia, serif'
      }
    },
    title: {
      style: {
        color: '#6f4327',
        fontWeight: '700',
        fontFamily: '"Cinzel", "Lora", Georgia, serif'
      }
    },
    subtitle: { style: { color: '#66503c' } },
    xAxis: { lineColor: '#b78f68', tickColor: '#b78f68', gridLineColor: 'rgba(111, 67, 39, 0.16)', labels: { style: { color: '#66503c' } } },
    yAxis: { lineColor: '#b78f68', tickColor: '#b78f68', gridLineColor: 'rgba(111, 67, 39, 0.16)', labels: { style: { color: '#66503c' } } },
    legend: { itemStyle: { color: '#4d3a2a' } },
    tooltip: { backgroundColor: '#f5e7d2', borderColor: '#b78f68', style: { color: '#2e2218' } },
    credits: { enabled: false }
  }

  var HIGHCHARTS_SCIFI1950S = {
    colors: ['#56f0ff', '#ff7a00', '#b7ff57', '#ff4f7a', '#9c88ff', '#79ffe1', '#ffd166'],
    chart: {
      backgroundColor: 'transparent',
      style: {
        color: '#d6f6ff',
        fontFamily: '"Orbitron", "Exo 2", "IBM Plex Sans", -apple-system, BlinkMacSystemFont, "Segoe UI", sans-serif'
      }
    },
    title: {
      style: {
        color: '#ffbe5c',
        fontWeight: '700',
        fontFamily: '"Orbitron", "Exo 2", "IBM Plex Sans", -apple-system, BlinkMacSystemFont, "Segoe UI", sans-serif'
      }
    },
    subtitle: { style: { color: '#8fdfff' } },
    xAxis: { lineColor: '#3eb6d6', tickColor: '#3eb6d6', gridLineColor: 'rgba(86, 240, 255, 0.2)', labels: { style: { color: '#8fdfff' } } },
    yAxis: { lineColor: '#3eb6d6', tickColor: '#3eb6d6', gridLineColor: 'rgba(86, 240, 255, 0.2)', labels: { style: { color: '#8fdfff' } } },
    legend: { itemStyle: { color: '#c1f0ff' } },
    tooltip: { backgroundColor: '#121e34', borderColor: '#56f0ff', style: { color: '#e9fbff' } },
    credits: { enabled: false }
  }

  var HIGHCHARTS_DRWHO = {
    colors: ['#4aa8ff', '#7b61ff', '#4de0ff', '#ffd166', '#ff6b9f', '#7ef29a', '#9bb8ff'],
    chart: {
      backgroundColor: 'transparent',
      style: {
        color: '#d9ecff',
        fontFamily: '"Orbitron", "Exo 2", "IBM Plex Sans", -apple-system, BlinkMacSystemFont, "Segoe UI", sans-serif'
      }
    },
    title: { style: { color: '#7ecbff', fontWeight: '700', fontFamily: '"Orbitron", "Exo 2", "IBM Plex Sans", -apple-system, BlinkMacSystemFont, "Segoe UI", sans-serif' } },
    subtitle: { style: { color: '#a8c8ff' } },
    xAxis: { lineColor: '#4aa8ff', tickColor: '#4aa8ff', gridLineColor: 'rgba(74, 168, 255, 0.2)', labels: { style: { color: '#a8c8ff' } } },
    yAxis: { lineColor: '#4aa8ff', tickColor: '#4aa8ff', gridLineColor: 'rgba(74, 168, 255, 0.2)', labels: { style: { color: '#a8c8ff' } } },
    legend: { itemStyle: { color: '#d9ecff' } },
    tooltip: { backgroundColor: '#101b33', borderColor: '#4aa8ff', style: { color: '#eef7ff' } },
    credits: { enabled: false }
  }

  var HIGHCHARTS_PLAYSCHOOL = {
    colors: ['#ff5e5e', '#ffd166', '#4ddf8e', '#4aa8ff', '#a66cff', '#ff8cc6', '#ff9f43'],
    chart: {
      backgroundColor: 'transparent',
      style: {
        color: '#28385a',
        fontFamily: '"Baloo 2", "Nunito", "Manrope", -apple-system, BlinkMacSystemFont, "Segoe UI", sans-serif'
      }
    },
    title: { style: { color: '#2f5fb8', fontWeight: '700', fontFamily: '"Baloo 2", "Nunito", "Manrope", sans-serif' } },
    subtitle: { style: { color: '#4c6492' } },
    xAxis: { lineColor: '#8eb8ff', tickColor: '#8eb8ff', gridLineColor: 'rgba(74, 168, 255, 0.16)', labels: { style: { color: '#4c6492' } } },
    yAxis: { lineColor: '#8eb8ff', tickColor: '#8eb8ff', gridLineColor: 'rgba(74, 168, 255, 0.16)', labels: { style: { color: '#4c6492' } } },
    legend: { itemStyle: { color: '#364a74' } },
    tooltip: { backgroundColor: '#f5f9ff', borderColor: '#8eb8ff', style: { color: '#28385a' } },
    credits: { enabled: false }
  }

  var HIGHCHARTS_TRON = {
    colors: ['#00f0ff', '#ff2fd0', '#00ffa8', '#7cf8ff', '#ff8a00', '#00b8ff', '#ffa6f0'],
    chart: { backgroundColor: 'transparent', style: { color: '#e8fdff', fontFamily: '"Orbitron", "Exo 2", "IBM Plex Sans", -apple-system, BlinkMacSystemFont, "Segoe UI", sans-serif' } },
    title: { style: { color: '#00f0ff', fontWeight: '700', fontFamily: '"Orbitron", "Exo 2", "IBM Plex Sans", -apple-system, BlinkMacSystemFont, "Segoe UI", sans-serif' } },
    subtitle: { style: { color: '#9efbff' } },
    xAxis: { lineColor: '#00f0ff', tickColor: '#00f0ff', gridLineColor: 'rgba(0, 240, 255, 0.35)', labels: { style: { color: '#a8f7ff' } } },
    yAxis: { lineColor: '#00f0ff', tickColor: '#00f0ff', gridLineColor: 'rgba(0, 240, 255, 0.35)', labels: { style: { color: '#a8f7ff' } } },
    legend: { itemStyle: { color: '#d4fcff' } },
    tooltip: { backgroundColor: '#050a12', borderColor: '#ff2fd0', style: { color: '#e8fdff' } },
    credits: { enabled: false }
  }

  var HIGHCHARTS_STARWARSANH = {
    colors: ['#ffe066', '#e5b800', '#7ad8ff', '#ff8f66', '#9be564', '#f6c85f', '#8bb6ff'],
    chart: { backgroundColor: 'transparent', style: { color: '#f6e8b1', fontFamily: '"Cinzel", "Lora", Georgia, serif' } },
    title: { style: { color: '#ffe066', fontWeight: '700', fontFamily: '"Cinzel", "Lora", Georgia, serif' } },
    subtitle: { style: { color: '#d9c98b' } },
    xAxis: { lineColor: '#d2b34a', tickColor: '#d2b34a', gridLineColor: 'rgba(255, 224, 102, 0.18)', labels: { style: { color: '#d9c98b' } } },
    yAxis: { lineColor: '#d2b34a', tickColor: '#d2b34a', gridLineColor: 'rgba(255, 224, 102, 0.18)', labels: { style: { color: '#d9c98b' } } },
    legend: { itemStyle: { color: '#f2ddb0' } },
    tooltip: { backgroundColor: '#15120b', borderColor: '#d2b34a', style: { color: '#f6e8b1' } },
    credits: { enabled: false }
  }

  var HIGHCHARTS_ACADEMIC = {
    colors: ['#1f4e79', '#2a6f97', '#5a8fbd', '#7c9f5b', '#b58b41', '#8f3f4f', '#6f7f92'],
    chart: { backgroundColor: 'transparent', style: { color: '#223243', fontFamily: '"IBM Plex Sans", "Source Sans 3", -apple-system, BlinkMacSystemFont, "Segoe UI", sans-serif' } },
    title: { style: { color: '#1f4e79', fontWeight: '700', fontFamily: '"IBM Plex Sans", "Source Sans 3", -apple-system, BlinkMacSystemFont, "Segoe UI", sans-serif' } },
    subtitle: { style: { color: '#4b6072' } },
    xAxis: { lineColor: '#9db4c8', tickColor: '#9db4c8', gridLineColor: 'rgba(31, 78, 121, 0.14)', labels: { style: { color: '#4b6072' } } },
    yAxis: { lineColor: '#9db4c8', tickColor: '#9db4c8', gridLineColor: 'rgba(31, 78, 121, 0.14)', labels: { style: { color: '#4b6072' } } },
    legend: { itemStyle: { color: '#31495f' } },
    tooltip: { backgroundColor: '#f6f9fc', borderColor: '#9db4c8', style: { color: '#223243' } },
    credits: { enabled: false }
  }

  var HIGHCHARTS_MUTEDSCIENTIST = {
    colors: ['#496b7a', '#6a8a92', '#8fa7ab', '#7f9271', '#9b8b75', '#7d6f82', '#64857a'],
    chart: { backgroundColor: 'transparent', style: { color: '#2a3438', fontFamily: '"Source Sans 3", "IBM Plex Sans", -apple-system, BlinkMacSystemFont, "Segoe UI", sans-serif' } },
    title: { style: { color: '#496b7a', fontWeight: '700', fontFamily: '"Source Sans 3", "IBM Plex Sans", -apple-system, BlinkMacSystemFont, "Segoe UI", sans-serif' } },
    subtitle: { style: { color: '#5f6f75' } },
    xAxis: { lineColor: '#b7c4c9', tickColor: '#b7c4c9', gridLineColor: 'rgba(73, 107, 122, 0.14)', labels: { style: { color: '#5f6f75' } } },
    yAxis: { lineColor: '#b7c4c9', tickColor: '#b7c4c9', gridLineColor: 'rgba(73, 107, 122, 0.14)', labels: { style: { color: '#5f6f75' } } },
    legend: { itemStyle: { color: '#475a62' } },
    tooltip: { backgroundColor: '#f1f5f6', borderColor: '#b7c4c9', style: { color: '#2a3438' } },
    credits: { enabled: false }
  }

  var HIGHCHARTS_TERMINAL = {
    colors: ['#33ff66', '#99ff33', '#66ff99', '#ffbf00', '#66ffcc', '#99ff66', '#ccff99'],
    chart: { backgroundColor: 'transparent', style: { color: '#9cff9c', fontFamily: '"VT323", "IBM Plex Sans", monospace' } },
    title: { style: { color: '#33ff66', fontWeight: '700', fontFamily: '"VT323", "IBM Plex Sans", monospace' } },
    subtitle: { style: { color: '#7fe27f' } },
    xAxis: { lineColor: '#33ff66', tickColor: '#33ff66', gridLineColor: 'rgba(51, 255, 102, 0.2)', labels: { style: { color: '#7fe27f' } } },
    yAxis: { lineColor: '#33ff66', tickColor: '#33ff66', gridLineColor: 'rgba(51, 255, 102, 0.2)', labels: { style: { color: '#7fe27f' } } },
    legend: { itemStyle: { color: '#9cff9c' } },
    tooltip: { backgroundColor: '#0a140a', borderColor: '#33ff66', style: { color: '#c8ffc8' } },
    credits: { enabled: false }
  }

  var HIGHCHARTS_KIMIKO = {
    colors: ['#d9292f', '#0f2a44', '#f2f2f0', '#1d3b5c', '#b71f28', '#7b8a99', '#e35d63'],
    chart: { backgroundColor: 'transparent', style: { color: '#1f2b37', fontFamily: '"Noto Sans JP", "IBM Plex Sans", "Source Sans 3", sans-serif' } },
    title: { style: { color: '#0f2a44', fontWeight: '700', fontFamily: '"Noto Sans JP", "IBM Plex Sans", "Source Sans 3", sans-serif' } },
    subtitle: { style: { color: '#4c5f72' } },
    xAxis: { lineColor: '#8ea1b2', tickColor: '#8ea1b2', gridLineColor: 'rgba(15, 42, 68, 0.12)', labels: { style: { color: '#4c5f72' } } },
    yAxis: { lineColor: '#8ea1b2', tickColor: '#8ea1b2', gridLineColor: 'rgba(15, 42, 68, 0.12)', labels: { style: { color: '#4c5f72' } } },
    legend: { itemStyle: { color: '#2d4054' } },
    tooltip: { backgroundColor: '#f7f7f4', borderColor: '#8ea1b2', style: { color: '#1f2b37' } },
    credits: { enabled: false }
  }

  var HIGHCHARTS_LONDONCALLING = {
    colors: ['#e83f7a', '#80f79a', '#f2f0e8', '#a6adb4', '#ff7aa8', '#9dffb6', '#d2d6db'],
    chart: { backgroundColor: 'transparent', style: { color: '#f2f0e8', fontFamily: '"Barlow Condensed", "IBM Plex Sans", "Source Sans 3", sans-serif' } },
    title: { style: { color: '#e83f7a', fontWeight: '700', fontFamily: '"Archivo Black", "Barlow Condensed", "IBM Plex Sans", sans-serif' } },
    subtitle: { style: { color: '#cfd4da' } },
    xAxis: { lineColor: '#7f8994', tickColor: '#7f8994', gridLineColor: 'rgba(242, 240, 232, 0.14)', labels: { style: { color: '#cfd4da' } } },
    yAxis: { lineColor: '#7f8994', tickColor: '#7f8994', gridLineColor: 'rgba(242, 240, 232, 0.14)', labels: { style: { color: '#cfd4da' } } },
    legend: { itemStyle: { color: '#f2f0e8' } },
    tooltip: { backgroundColor: '#121315', borderColor: '#e83f7a', style: { color: '#f2f0e8' } },
    credits: { enabled: false }
  }

  var HIGHCHARTS_GRUMPYCHEMIST = {
    colors: ['#c9a227', '#6b9b7a', '#e8e6e0', '#8b929c', '#e4c86a', '#8fbc9a', '#5c6570'],
    chart: { backgroundColor: 'transparent', style: { color: '#e8e6e0', fontFamily: '"IBM Plex Sans", "Source Sans 3", sans-serif' } },
    title: { style: { color: '#d4b44a', fontWeight: '700', fontFamily: '"Libre Baskerville", "Lora", Georgia, serif' } },
    subtitle: { style: { color: '#a8b0ba' } },
    xAxis: { lineColor: '#5c6570', tickColor: '#5c6570', gridLineColor: 'rgba(148, 163, 184, 0.14)', labels: { style: { color: '#a8b0ba' } } },
    yAxis: { lineColor: '#5c6570', tickColor: '#5c6570', gridLineColor: 'rgba(148, 163, 184, 0.14)', labels: { style: { color: '#a8b0ba' } } },
    legend: { itemStyle: { color: '#e8e6e0' } },
    tooltip: { backgroundColor: '#1e2329', borderColor: '#c9a227', style: { color: '#e8e6e0' } },
    credits: { enabled: false }
  }

  var HIGHCHARTS_OXFORD = {
    colors: ['#002147', '#8b6914', '#1c1b18', '#003875', '#4a6fa5', '#c4b8a8', '#2d4a6f'],
    chart: { backgroundColor: 'transparent', style: { color: '#1c1b18', fontFamily: '"EB Garamond", "Cormorant Garamond", Georgia, serif' } },
    title: { style: { color: '#002147', fontWeight: '700', fontFamily: '"Cinzel", "Libre Baskerville", Georgia, serif' } },
    subtitle: { style: { color: '#4a4a42' } },
    xAxis: { lineColor: '#002147', tickColor: '#002147', gridLineColor: 'rgba(0, 33, 71, 0.12)', labels: { style: { color: '#3d3d36' } } },
    yAxis: { lineColor: '#002147', tickColor: '#002147', gridLineColor: 'rgba(0, 33, 71, 0.12)', labels: { style: { color: '#3d3d36' } } },
    legend: { itemStyle: { color: '#1c1b18' } },
    tooltip: { backgroundColor: '#faf6ef', borderColor: '#002147', style: { color: '#1c1b18' } },
    credits: { enabled: false }
  }

  var HIGHCHARTS_CAMBRIDGE = {
    colors: ['#0ea5e9', '#0284c7', '#0369a1', '#38bdf8', '#7dd3fc', '#bae6fd', '#0c4a6e'],
    chart: { backgroundColor: 'transparent', style: { color: '#0c4a6e', fontFamily: '"Manrope", "IBM Plex Sans", "Source Sans 3", sans-serif' } },
    title: { style: { color: '#0369a1', fontWeight: '700', fontFamily: '"Manrope", "IBM Plex Sans", sans-serif' } },
    subtitle: { style: { color: '#64748b' } },
    xAxis: { lineColor: '#7dd3fc', tickColor: '#7dd3fc', gridLineColor: 'rgba(14, 165, 233, 0.15)', labels: { style: { color: '#475569' } } },
    yAxis: { lineColor: '#7dd3fc', tickColor: '#7dd3fc', gridLineColor: 'rgba(14, 165, 233, 0.15)', labels: { style: { color: '#475569' } } },
    legend: { itemStyle: { color: '#0c4a6e' } },
    tooltip: { backgroundColor: '#ffffff', borderColor: '#38bdf8', style: { color: '#0c4a6e' } },
    credits: { enabled: false }
  }

  function getThemePalette (theme) {
    if (theme === 'dark') {
      return {
        plotBg: '#0f172a',
        heatmapMin: '#0f172a',
        heatmapMax: '#5eead4',
        heatmapNull: '#1f2937',
        heatmapBorder: 'rgba(148, 163, 184, 0.22)',
        metaColor: '#0f172a'
      }
    }
    if (theme === 'hacienda') {
      return {
        plotBg: '#1a0d26',
        heatmapMin: '#1a0d26',
        heatmapMax: '#ff3b8d',
        heatmapNull: '#2f1841',
        heatmapBorder: 'rgba(255, 127, 17, 0.35)',
        metaColor: '#1a0d26'
      }
    }
    if (theme === 'sf1950s') {
      return {
        plotBg: '#f4ecde',
        heatmapMin: '#f4ecde',
        heatmapMax: '#4f7a78',
        heatmapNull: '#d9cfbf',
        heatmapBorder: '#c8b79f',
        metaColor: '#f4ecde'
      }
    }
    if (theme === 'paris1940s') {
      return {
        plotBg: '#efe1cb',
        heatmapMin: '#efe1cb',
        heatmapMax: '#7a4f2f',
        heatmapNull: '#d4c4ab',
        heatmapBorder: '#bba78c',
        metaColor: '#efe1cb'
      }
    }
    if (theme === 'hailmary') {
      return {
        plotBg: '#0b1322',
        heatmapMin: '#0b1322',
        heatmapMax: '#20d8ff',
        heatmapNull: '#14243a',
        heatmapBorder: 'rgba(125, 249, 255, 0.25)',
        metaColor: '#0b1322'
      }
    }
    if (theme === 'stoneage') {
      return {
        plotBg: '#efe0c4',
        heatmapMin: '#efe0c4',
        heatmapMax: '#b3532e',
        heatmapNull: '#d9c3a0',
        heatmapBorder: '#c79c64',
        metaColor: '#efe0c4'
      }
    }
    if (theme === 'disney') {
      return {
        plotBg: '#eaf5ff',
        heatmapMin: '#eaf5ff',
        heatmapMax: '#1e6fd9',
        heatmapNull: '#cfe6fb',
        heatmapBorder: '#9fc8f8',
        metaColor: '#eaf5ff'
      }
    }
    if (theme === 'pixar') {
      return {
        plotBg: '#edf7ff',
        heatmapMin: '#edf7ff',
        heatmapMax: '#0b5fff',
        heatmapNull: '#d8ecff',
        heatmapBorder: '#8bc5ff',
        metaColor: '#edf7ff'
      }
    }
    if (theme === 'telextext') {
      return {
        plotBg: '#05101c',
        heatmapMin: '#05101c',
        heatmapMax: '#00ffff',
        heatmapNull: '#0c1e2f',
        heatmapBorder: 'rgba(0, 255, 255, 0.25)',
        metaColor: '#05101c'
      }
    }
    if (theme === 'flowerpower1960s') {
      return {
        plotBg: '#fff6d6',
        heatmapMin: '#fff6d6',
        heatmapMax: '#ff4fa1',
        heatmapNull: '#f8e1ff',
        heatmapBorder: '#d77bff',
        metaColor: '#fff6d6'
      }
    }
    if (theme === 'marblegold') {
      return {
        plotBg: '#f6f6f3',
        heatmapMin: '#f6f6f3',
        heatmapMax: '#c8a24a',
        heatmapNull: '#e7e6e1',
        heatmapBorder: '#c8b27a',
        metaColor: '#f6f6f3'
      }
    }
    if (theme === 'rusticwood') {
      return {
        plotBg: '#f3e1c8',
        heatmapMin: '#f3e1c8',
        heatmapMax: '#8c5a3c',
        heatmapNull: '#e4ccb0',
        heatmapBorder: '#b78f68',
        metaColor: '#f3e1c8'
      }
    }
    if (theme === 'scifi1950s') {
      return {
        plotBg: '#121e34',
        heatmapMin: '#121e34',
        heatmapMax: '#56f0ff',
        heatmapNull: '#1d2b45',
        heatmapBorder: 'rgba(86, 240, 255, 0.25)',
        metaColor: '#121e34'
      }
    }
    if (theme === 'drwho') {
      return {
        plotBg: '#101b33',
        heatmapMin: '#101b33',
        heatmapMax: '#4aa8ff',
        heatmapNull: '#1b2742',
        heatmapBorder: 'rgba(74, 168, 255, 0.25)',
        metaColor: '#101b33'
      }
    }
    if (theme === 'playschool') {
      return {
        plotBg: '#f2f8ff',
        heatmapMin: '#f2f8ff',
        heatmapMax: '#4aa8ff',
        heatmapNull: '#dfeeff',
        heatmapBorder: '#8eb8ff',
        metaColor: '#f2f8ff'
      }
    }
    if (theme === 'tron') {
      return {
        plotBg: '#050a12',
        heatmapMin: '#050a12',
        heatmapMax: '#00f0ff',
        heatmapNull: '#0a1828',
        heatmapBorder: 'rgba(0, 240, 255, 0.4)',
        metaColor: '#050a12'
      }
    }
    if (theme === 'starwarsanh') {
      return {
        plotBg: '#15120b',
        heatmapMin: '#15120b',
        heatmapMax: '#ffe066',
        heatmapNull: '#2a2417',
        heatmapBorder: 'rgba(210, 179, 74, 0.28)',
        metaColor: '#15120b'
      }
    }
    if (theme === 'academic') {
      return {
        plotBg: '#f6f9fc',
        heatmapMin: '#f6f9fc',
        heatmapMax: '#1f4e79',
        heatmapNull: '#e7eef5',
        heatmapBorder: '#9db4c8',
        metaColor: '#f6f9fc'
      }
    }
    if (theme === 'mutedscientist') {
      return {
        plotBg: '#f1f5f6',
        heatmapMin: '#f1f5f6',
        heatmapMax: '#496b7a',
        heatmapNull: '#e2e9ec',
        heatmapBorder: '#b7c4c9',
        metaColor: '#f1f5f6'
      }
    }
    if (theme === 'terminal') {
      return {
        plotBg: '#0a140a',
        heatmapMin: '#0a140a',
        heatmapMax: '#33ff66',
        heatmapNull: '#132613',
        heatmapBorder: 'rgba(51, 255, 102, 0.25)',
        metaColor: '#0a140a'
      }
    }
    if (theme === 'kimiko') {
      return {
        plotBg: '#f7f7f4',
        heatmapMin: '#f7f7f4',
        heatmapMax: '#d9292f',
        heatmapNull: '#e6e7e5',
        heatmapBorder: '#8ea1b2',
        metaColor: '#f7f7f4'
      }
    }
    if (theme === 'londoncalling') {
      return {
        plotBg: '#121315',
        heatmapMin: '#121315',
        heatmapMax: '#e83f7a',
        heatmapNull: '#1d1f23',
        heatmapBorder: 'rgba(242, 240, 232, 0.22)',
        metaColor: '#121315'
      }
    }
    if (theme === 'grumpychemist') {
      return {
        plotBg: '#1e2329',
        heatmapMin: '#1e2329',
        heatmapMax: '#c9a227',
        heatmapNull: '#2a3038',
        heatmapBorder: 'rgba(201, 162, 39, 0.25)',
        metaColor: '#1e2329'
      }
    }
    if (theme === 'oxford') {
      return {
        plotBg: '#f4f0e6',
        heatmapMin: '#f4f0e6',
        heatmapMax: '#002147',
        heatmapNull: '#e0d8cc',
        heatmapBorder: 'rgba(0, 33, 71, 0.35)',
        metaColor: '#f2ebe0'
      }
    }
    if (theme === 'cambridge') {
      return {
        plotBg: '#f0f9ff',
        heatmapMin: '#f0f9ff',
        heatmapMax: '#0284c7',
        heatmapNull: '#e0f2fe',
        heatmapBorder: 'rgba(56, 189, 248, 0.35)',
        metaColor: '#e0f2fe'
      }
    }
    return {
      plotBg: '#ffffff',
      heatmapMin: '#ffffff',
      heatmapMax: '#004B44',
      heatmapNull: '#EFEFEF',
      heatmapBorder: '#e2e8f0',
      metaColor: '#ffffff'
    }
  }

  function getStoredTheme () {
    try {
      var t = window.localStorage.getItem(STORAGE_KEY)
      if (VALID_THEMES.indexOf(t) !== -1) return t
    } catch (e) {}
    return null
  }

  /**
   * Charts created before a theme switch keep old plotBackgroundColor; update in place.
   */
  function reflowExistingCharts (theme) {
    if (typeof window.Highcharts === 'undefined') return
    var palette = getThemePalette(theme)
    var charts = window.Highcharts.charts
    if (!charts || !charts.length) return
    for (var i = 0; i < charts.length; i++) {
      var chart = charts[i]
      if (!chart || typeof chart.update !== 'function') continue
      var hasHeatmap = false
      try {
        hasHeatmap = !!(chart.series && chart.series.some(function (s) { return s && s.type === 'heatmap' }))
      } catch (e) {}
      try {
        var updateOptions = {
          chart: {
            backgroundColor: 'transparent',
            plotBackgroundColor: palette.plotBg
          }
        }
        if (hasHeatmap) {
          updateOptions.colorAxis = {
            minColor: palette.heatmapMin,
            maxColor: palette.heatmapMax
          }
          updateOptions.plotOptions = {
            heatmap: {
              nullColor: palette.heatmapNull,
              borderColor: palette.heatmapBorder,
              borderWidth: 1
            }
          }
        }
        chart.update(updateOptions, true)
      } catch (e) {}
    }
  }

  function applyHighcharts (theme) {
    if (typeof window.Highcharts === 'undefined') return
    if (theme === 'dark') {
      window.Highcharts.setOptions(HIGHCHARTS_DARK)
    } else if (theme === 'hacienda') {
      window.Highcharts.setOptions(HIGHCHARTS_HACIENDA)
    } else if (theme === 'sf1950s') {
      window.Highcharts.setOptions(HIGHCHARTS_SF1950S)
    } else if (theme === 'paris1940s') {
      window.Highcharts.setOptions(HIGHCHARTS_PARIS1940S)
    } else if (theme === 'hailmary') {
      window.Highcharts.setOptions(HIGHCHARTS_HAILMARY)
    } else if (theme === 'stoneage') {
      window.Highcharts.setOptions(HIGHCHARTS_STONEAGE)
    } else if (theme === 'disney') {
      window.Highcharts.setOptions(HIGHCHARTS_DISNEY)
    } else if (theme === 'pixar') {
      window.Highcharts.setOptions(HIGHCHARTS_PIXAR)
    } else if (theme === 'telextext') {
      window.Highcharts.setOptions(HIGHCHARTS_TELEXTEXT)
    } else if (theme === 'flowerpower1960s') {
      window.Highcharts.setOptions(HIGHCHARTS_FLOWERPOWER1960S)
    } else if (theme === 'marblegold') {
      window.Highcharts.setOptions(HIGHCHARTS_MARBLEGOLD)
    } else if (theme === 'rusticwood') {
      window.Highcharts.setOptions(HIGHCHARTS_RUSTICWOOD)
    } else if (theme === 'scifi1950s') {
      window.Highcharts.setOptions(HIGHCHARTS_SCIFI1950S)
    } else if (theme === 'drwho') {
      window.Highcharts.setOptions(HIGHCHARTS_DRWHO)
    } else if (theme === 'playschool') {
      window.Highcharts.setOptions(HIGHCHARTS_PLAYSCHOOL)
    } else if (theme === 'tron') {
      window.Highcharts.setOptions(HIGHCHARTS_TRON)
    } else if (theme === 'starwarsanh') {
      window.Highcharts.setOptions(HIGHCHARTS_STARWARSANH)
    } else if (theme === 'academic') {
      window.Highcharts.setOptions(HIGHCHARTS_ACADEMIC)
    } else if (theme === 'mutedscientist') {
      window.Highcharts.setOptions(HIGHCHARTS_MUTEDSCIENTIST)
    } else if (theme === 'terminal') {
      window.Highcharts.setOptions(HIGHCHARTS_TERMINAL)
    } else if (theme === 'kimiko') {
      window.Highcharts.setOptions(HIGHCHARTS_KIMIKO)
    } else if (theme === 'londoncalling') {
      window.Highcharts.setOptions(HIGHCHARTS_LONDONCALLING)
    } else if (theme === 'grumpychemist') {
      window.Highcharts.setOptions(HIGHCHARTS_GRUMPYCHEMIST)
    } else if (theme === 'oxford') {
      window.Highcharts.setOptions(HIGHCHARTS_OXFORD)
    } else if (theme === 'cambridge') {
      window.Highcharts.setOptions(HIGHCHARTS_CAMBRIDGE)
    } else if (window.__minotourHcLightTheme) {
      window.Highcharts.setOptions(window.__minotourHcLightTheme)
    }
    reflowExistingCharts(theme)
    try {
      window.dispatchEvent(new CustomEvent('minotour:themechange', { detail: { theme: theme } }))
    } catch (e) {}
  }

  function setMetaThemeColor (theme) {
    var meta = document.querySelector('meta[name="theme-color"]')
    if (!meta) return
    meta.setAttribute('content', getThemePalette(theme).metaColor)
  }

  function syncDarkUiClass (theme) {
    if (DARK_UI_THEMES[theme]) document.documentElement.classList.add('dark')
    else document.documentElement.classList.remove('dark')
  }

  function setTheme (theme) {
    if (VALID_THEMES.indexOf(theme) === -1) theme = 'light'
    document.documentElement.setAttribute('data-theme', theme)
    syncDarkUiClass(theme)
    try {
      window.localStorage.setItem(STORAGE_KEY, theme)
    } catch (e) {}
    applyHighcharts(theme)
    setMetaThemeColor(theme)
    syncToggleUi()
  }

  function syncToggleUi () {
    var root = document.documentElement
    var theme = root.getAttribute('data-theme') || 'light'
    var isDark = theme === 'dark'
    var select = document.getElementById('minotour-theme-select')
    if (select && select.value !== theme) select.value = theme
    var btn = document.getElementById('minotour-theme-toggle')
    if (btn) {
      btn.setAttribute('aria-checked', isDark ? 'true' : 'false')
      var on = btn.querySelector('[data-minotour-theme-state="on"]')
      var off = btn.querySelector('[data-minotour-theme-state="off"]')
      if (on) on.classList.toggle('tw-hidden', !isDark)
      if (off) off.classList.toggle('tw-hidden', isDark)
    }
  }

  function init () {
    var stored = getStoredTheme()
    if (stored) {
      document.documentElement.setAttribute('data-theme', stored)
    }
    var t = document.documentElement.getAttribute('data-theme') || 'light'
    if (VALID_THEMES.indexOf(t) === -1) t = 'light'
    syncDarkUiClass(t)
    applyHighcharts(t)
    setMetaThemeColor(t)

    var select = document.getElementById('minotour-theme-select')
    if (select) {
      select.addEventListener('change', function (e) {
        setTheme(e.target.value)
      })
    }

    var btn = document.getElementById('minotour-theme-toggle')
    if (btn) {
      btn.addEventListener('click', function () {
        var cur = document.documentElement.getAttribute('data-theme') || 'light'
        setTheme(cur === 'dark' ? 'light' : 'dark')
      })
    }
    syncToggleUi()
  }

  if (document.readyState === 'loading') {
    document.addEventListener('DOMContentLoaded', init)
  } else {
    init()
  }

  window.MinoTourTheme = {
    STORAGE_KEY: STORAGE_KEY,
    setTheme: setTheme,
    getTheme: function () {
      return document.documentElement.getAttribute('data-theme') || 'light'
    },
    reflowCharts: function () {
      reflowExistingCharts(document.documentElement.getAttribute('data-theme') || 'light')
    }
  }
})(window)
