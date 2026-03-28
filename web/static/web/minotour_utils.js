Highcharts.setOptions({
  plotOptions: {
    series: {
      animation: false,
      turboThreshold: 0
    }
  }
})

/**
 * Get the selected flowcell tab
 */
function getSelectedTab () {
  return sessionStorage.getItem(`flowcellTab`)
}

/**
 * Set the selected flowcell Tab
 * @param tab {string} The tab that a user has just switched to
 */
function setSelectedTab (tab) {
  sessionStorage.setItem(`flowcellTab`, tab)
}

/**
 * Set the Flowcell ID in the session storage for whatever flowcell we are inspecting.
 * @param flowcell_id
 */
function setSelectedFlowcell (flowcellId) {
  sessionStorage.setItem(`flowcellId`, flowcellId)
}

/**
 * Return the selected Flowcell from the hidden input.
 */
function getSelectedFlowcell () {
  return sessionStorage.getItem(`flowcellId`)
}

/**
 * Return the barcode name for given tab from session storage.
 * @param tab {string} Tab name
 * @returns {string} Barcode name for this tab.
 */
function getSelectedBarcode (tab) {
  return sessionStorage.getItem(`${tab}-barcode`)
}

/**
 * Set the barcode of the tab the user is on in the session storage.
 * @param barcode {string} Barcode name
 * @param tab {string} Tab name
 */
function setSelectedBarcode (barcode, tab) {
  sessionStorage.setItem(`${tab}-barcode`, barcode)
}

/**
 * Remove old series from a given chart
 * @param chart {object}
 */
function clearChartData (chart) {
  while (chart.series.length) {
    chart.series[0].remove(false)
  }
}

/**
 * Check whether the data we are about to add to a chart series is identical to the data already in that series.
 * @param newChartData {array} The data we plan to insert.
 * @param oldChartData {array} The data that exists in the series. Accessed via series.options.data.
 * @returns {bool} Whether the arrays are identical or not.
 */
function checkHighChartsDataIsIdentical (newChartData, oldChartData) {
  let flattenedNewData = []
  let flattenedOldData = []
  let identical = false
  flattenedNewData = [].concat(...newChartData)
  // flatten chart data so we can compare it to the new data, can't compare nested arrays
  flattenedOldData = [].concat(...oldChartData)

  // compare old data with new data
  identical = flattenedNewData.length === flattenedOldData.length && flattenedNewData.every((value, index) => value === flattenedOldData[index])
  return identical
}

/**
 * Get a cookie value.
 * @param name {string} The name of the cookie to be parsed.
 * @returns {string} The value of the given cookie.
 */
function getCookie (name) {
  let cookie
  let cookies
  let cookieValue = null
  if (document.cookie && document.cookie !== ``) {
    cookies = document.cookie.split(`;`)
    for (var i = 0; i < cookies.length; i++) {
      cookie = jQuery.trim(cookies[i])
      // Does this cookie string begin with the name we want?
      if (cookie.substring(0, name.length + 1) === (`${name}=`)) {
        cookieValue = decodeURIComponent(cookie.substring(name.length + 1))
        break
      }
    }
  }
  return cookieValue
}

/**
 * Set the active navbar item on the top main navbar
 * @param {number} item_index The index of the navbar element in the
 * @function {function} Set the active item on the top navbar.
 */
function setActiveNavbarItem (item_index) {
  var nav_bar = document.querySelectorAll(`.navbar li a`)

  nav_bar.forEach(function (element, index) {
    element.classList.remove(`active`)
    if (index == item_index) element.classList.add(`active`)
  })
}

/**
 * @function makeColumnChart
 * @param {string} divId ID of the div we want to draw the charts in
 * @param {string} chartTitle
 * @param {string} yAxisTitle
 * Creates an empty highCharts column chart at the provided divName,
 * with divName being the Id if the div we want to draw in.
 */
function makeColumnChart (divId, chartTitle, yAxisTitle) {
  return Highcharts.chart(divId, {
    chart: {
      type: `column`,
      animation: Highcharts.svg, // don"t animate in old IE
      marginRight: 10,
      zoomType: `x`
    },
    boost: {
      seriesThreshold: 10
    },
    title: {
      text: chartTitle
    },
    yAxis: {
      title: {
        text: yAxisTitle
      },
      plotLines: [{
        value: 0,
        width: 1,
        color: `#cbd5e1`
      }]
    },
    lang: {
      noData: `Looking for data probably.`
    },
    legend: {
      enabled: true
    },
    exporting: {
      enabled: true
    },
    series: []
  })
}

/**
 * Function abstracting creation of a HighCharts Spline chart.
 * @param divId {string} ID of div to initialise chart inside.
 * @param chartTitle {string} Title of chart.
 * @param yAxisTitle {string} Title of Y axis.
 * @returns {*} HighCharts Api object
 */
function makeSplineChart (divId, chartTitle, yAxisTitle) {
  const chart = Highcharts.chart(divId, {
    chart: {
      type: `spline`,
      marginRight: 10,
      animation: false,
      zoomType: `x`
    },
    boost: {
      enabled: false,
      useGPUTranslations: true
    },
    title: {
      text: chartTitle
    },
    xAxis: {
      type: `datetime`,
      tickPixelInterval: 150
    },
    yAxis: {
      title: {
        text: yAxisTitle
      },
      plotLines: [{
        value: 0,
        width: 1,
        color: `#cbd5e1`
      }]
    },
    legend: {
      enabled: true
    },
    exporting: {
      enabled: true,
      sourceWidth: 1200,
      sourceHeight: 400
    }
  })

  return chart
}

/**
 * Function abstracting creation of a HighCharts Spline chart without date time slider.
 * @param divName {string} ID of div to initialise chart inside.
 * @param chartTitle {string} Title of chart.
 * @param yAxisTitle {string} Title of Y axis.
 * @param xAxisTitle {string} Title of the X axis
 * @return {*} Highcharts API object
 */
function makeSplineChartNonDatetime (divName, chartTitle, yAxisTitle, xAxisTitle) {
  return Highcharts.chart(divName, {
    chart: {
      type: `spline`,
      marginRight: 10,
      animation: false,
      zoomType: `x`
    },
    boost: {
      useGPUTranslations: true
    },
    title: {
      text: chartTitle
    },
    xAxis: {
      title: {
        text: xAxisTitle
      }
    },
    yAxis: {
      title: {
        text: yAxisTitle
      },
      plotLines: [{
        value: 0,
        width: 1,
        color: `#cbd5e1`
      }]
    },
    legend: {
      enabled: true
    },
    exporting: {
      enabled: true
    }
  })
}

/**
 * TODO docstrin
 * @param divName
 * @param chartTitle
 * @param yAxisTitle
 */
function makeBoxPlot (divName, chartTitle, yAxisTitle) {
  return Highcharts.chart(divName, {
    chart: {
      type: `boxplot`
    },
    title: {
      text: chartTitle
    },
    legend: {
      enabled: true
    },
    xAxis: {
      type: `category`
    },
    yAxis: {
      title: {
        text: yAxisTitle
      },
      min: 0
    }
  })
}

/**
 * Read --mtc-* tokens from .flowcell-detail-page (obsidian-theme-compat.css) for Live Event charts.
 */
function getLiveChartsThemeFromDom () {
  var page = document.querySelector(`.flowcell-detail-page`)
  var cs = page ? getComputedStyle(page) : null
  function pick (name, fallback) {
    if (!cs) return fallback
    var v = cs.getPropertyValue(name).trim()
    return (v && v.length) ? v : fallback
  }
  return {
    surface: pick(`--mtc-surface`, `#ffffff`),
    surfaceAlt: pick(`--mtc-surface-alt`, `#f8fafc`),
    header: pick(`--mtc-header`, `#f1f5f9`),
    text: pick(`--mtc-text`, `#0f172a`),
    muted: pick(`--mtc-muted`, `#64748b`),
    border: pick(`--mtc-border`, `#cbd5e1`),
    pass: pick(`--mtc-pass`, `rgba(16, 185, 129, 0.12)`),
    accent: pick(`--obsidian`, `#004B44`)
  }
}

function mergeLiveStockChartTheme () {
  var t = getLiveChartsThemeFromDom()
  return {
    chart: {
      backgroundColor: t.surface,
      plotBackgroundColor: t.surface,
      plotBorderWidth: 1,
      plotBorderColor: t.border
    },
    tooltip: {
      backgroundColor: t.surfaceAlt,
      borderColor: t.border,
      style: {
        color: t.text
      }
    },
    rangeSelector: {
      buttonTheme: {
        fill: t.surfaceAlt,
        stroke: t.border,
        'stroke-width': 1,
        r: 2,
        style: {
          color: t.text
        },
        states: {
          hover: {
            fill: t.header
          },
          select: {
            fill: t.accent,
            style: {
              color: `#ffffff`
            }
          }
        }
      },
      inputBoxBorderColor: t.border,
      inputStyle: {
        color: t.text,
        fontWeight: `600`
      },
      labelStyle: {
        color: t.muted,
        fontWeight: `600`
      }
    },
    navigator: {
      outlineColor: t.border,
      maskFill: `rgba(0, 75, 68, 0.1)`,
      series: {
        color: t.accent,
        lineColor: t.accent
      }
    },
    legend: {
      itemStyle: { color: t.muted },
      itemHoverStyle: { color: t.text }
    },
    title: {
      style: {
        color: t.text,
        fontSize: `16px`,
        fontWeight: `700`
      }
    },
    xAxis: {
      lineColor: t.border,
      tickColor: t.border,
      labels: { style: { color: t.muted } },
      gridLineColor: t.border
    },
    yAxis: {
      lineColor: t.border,
      tickColor: t.border,
      gridLineColor: t.border,
      labels: { style: { color: t.muted } },
      title: { style: { color: t.muted } }
    }
  }
}

function mergeLiveColumnChartTheme () {
  var t = getLiveChartsThemeFromDom()
  return {
    chart: {
      backgroundColor: t.surface,
      plotBackgroundColor: t.surface,
      plotBorderWidth: 1,
      plotBorderColor: t.border
    },
    tooltip: {
      backgroundColor: t.surfaceAlt,
      borderColor: t.border,
      style: { color: t.text }
    },
    title: {
      style: {
        color: t.text,
        fontSize: `16px`,
        fontWeight: `700`
      }
    },
    legend: {
      itemStyle: { color: t.muted }
    },
    xAxis: {
      lineColor: t.border,
      labels: { style: { color: t.muted } },
      gridLineColor: t.border
    },
    yAxis: {
      lineColor: t.border,
      gridLineColor: t.border,
      labels: { style: { color: t.muted } },
      title: { style: { color: t.muted } }
    },
    plotOptions: {
      column: {
        borderWidth: 0,
        color: t.accent
      }
    }
  }
}

/**
 *
 * @param divId
 * @param chartTitle
 * @param yAxisTitle
 * @return {*}
 */
function makeLiveHistogram (divId, chartTitle, yAxisTitle) {
  var th = mergeLiveColumnChartTheme()
  return Highcharts.chart(divId, Highcharts.merge(true, th, {
    chart: {
      type: `column`,
      animation: false,
      zoomType: `x`,
      height: 400
    },
    title: {
      text: chartTitle
    },
    xAxis: {
      categories: []
    },
    yAxis: {
      title: {
        text: `Total Event Length`
      }
    },
    credits: {
      enabled: false
    },
    series: [{
      name: `Read Histogram`
    }]
  }))
}

/**
 *
 * @param divName {string} The name of the Div that we are drawing the chart in
 * @param chartTitle {string} The Title of the area chart
 * @param yAxisTitle {string} The y axis title of the area chart
 */
function makeAreaPlot (divName, chartTitle, yAxisTitle) {
  var th = mergeLiveStockChartTheme()
  return Highcharts.stockChart(divName, Highcharts.merge(true, th, {
    chart: {
      type: `area`,
      height: 350,
      marginRight: 200
    },
    boost: {
      useGPUTranslations: true
    },
    title: {
      text: chartTitle
    },
    xAxis: {
      range: 1 * 360 * 1000
    },
    rangeSelector: {
      enabled: true,
      selected: 2
    },
    yAxis: {
      endOnTick: false,
      title: {
        text: `Channel Classifications`
      }
    },
    legend: {
      enabled: true,
      align: `right`,
      verticalAlign: `top`,
      layout: `vertical`,
      x: 10
    },
    plotOptions: {
      area: {
        stacking: `percent`
      },
      series: {
        showInNavigator: true,
        dataLabels: {
          enabled: false,
          formatter: function () {
            return this.y
          }
        }
      }
    },
    credits: {
      enabled: false
    },
    series: []
  }))
}

// Todo can these spline chart creation functions be merged??
/**
 * Crate a spline chart with date range slider
 * @param divId {string} The id of the div that we are appending to
 * @param chartTitle {string} The chart title to be displayed
 * @param yAxisTitle {string} The chart yAxis title to be displayed
 * @return {*} HighCharts api object
 */
function makeLiveChart (divId, chartTitle, yAxisTitle) {
  var t = getLiveChartsThemeFromDom()
  var th = mergeLiveStockChartTheme()
  /* Fixed height: charts are constructed while the live tab is still .hidden, so container
   * offsetHeight is 0 and Highcharts would otherwise measure a collapsed box. */
  return Highcharts.stockChart(divId, Highcharts.merge(true, th, {
    chart: {
      type: `spline`,
      zoomType: `x`,
      height: 420
    },
    boost: {
      useGPUTranslations: true
    },
    rangeSelector: {
      enabled: true,
      buttons: [{
        type: `minute`,
        count: 1,
        text: `1min`
      }, {
        type: `minute`,
        count: 5,
        text: `5min`
      }, {
        type: `minute`,
        count: 30,
        text: `1/2hr`
      }, {
        type: `minute`,
        count: 60,
        text: `1hr`
      }, {
        type: `day`,
        count: 0.5,
        text: `12hrs`
      }, {
        type: `day`,
        count: 1,
        text: `1day`
      }, {
        type: `all`,
        text: `All`
      }]
    },
    title: {
      text: chartTitle
    },
    xAxis: {
      type: `datetime`,
      tickPixelInterval: 150
    },
    yAxis: {
      title: {
        text: yAxisTitle
      },
      plotLines: [{
        value: 0,
        width: 1,
        color: t.border
      }]
    },
    credits: {
      enabled: false
    },
    series: []
  }))
};

/**
 * Create a HighCharts heatmap chart. These are used on the basecalled data tab.
 * @param divName {string} The id of the div to render chart to.
 * @param chartTitle {string} The title of the chart.
 * @return {*}
 */
function makeHeatmapChart (divName, chartTitle) {
  var theme = document.documentElement.getAttribute('data-theme') || 'light'
  var isDark = theme === 'dark'
  var isHacienda = theme === 'hacienda'
  var isSF1950s = theme === 'sf1950s'
  var isParis1940s = theme === 'paris1940s'
  var isHailMary = theme === 'hailmary'
  var isStoneage = theme === 'stoneage'
  var isDisney = theme === 'disney'
  var isPixar = theme === 'pixar'
  var isTelextext = theme === 'telextext'
  var isFlowerPower1960s = theme === 'flowerpower1960s'
  var isMarbleGold = theme === 'marblegold'
  var isRusticWood = theme === 'rusticwood'
  var isScifi1950s = theme === 'scifi1950s'
  var isDrWho = theme === 'drwho'
  var isPlaySchool = theme === 'playschool'
  var isTron = theme === 'tron'
  var isStarWarsANH = theme === 'starwarsanh'
  var isAcademic = theme === 'academic'
  var isMutedScientist = theme === 'mutedscientist'
  var isTerminal = theme === 'terminal'
  var isKimiko = theme === 'kimiko'
  var isLondonCalling = theme === 'londoncalling'
  var isGrumpyChemist = theme === 'grumpychemist'
  var isOxford = theme === 'oxford'
  var isCambridge = theme === 'cambridge'
  var plotBg = isDark ? '#0f172a' : (isHacienda ? '#1a0d26' : (isSF1950s ? '#f4ecde' : (isParis1940s ? '#efe1cb' : (isHailMary ? '#0b1322' : (isStoneage ? '#efe0c4' : (isDisney ? '#eaf5ff' : (isPixar ? '#edf7ff' : (isTelextext ? '#05101c' : (isFlowerPower1960s ? '#fff6d6' : (isMarbleGold ? '#f6f6f3' : (isRusticWood ? '#f3e1c8' : '#ffffff')))))))))))
  var minColor = isDark ? '#0f172a' : (isHacienda ? '#1a0d26' : (isSF1950s ? '#f4ecde' : (isParis1940s ? '#efe1cb' : (isHailMary ? '#0b1322' : (isStoneage ? '#efe0c4' : (isDisney ? '#eaf5ff' : (isPixar ? '#edf7ff' : (isTelextext ? '#05101c' : (isFlowerPower1960s ? '#fff6d6' : (isMarbleGold ? '#f6f6f3' : (isRusticWood ? '#f3e1c8' : '#FFFFFF')))))))))))
  var maxColor = isDark ? '#5eead4' : (isHacienda ? '#ff3b8d' : (isSF1950s ? '#4f7a78' : (isParis1940s ? '#7a4f2f' : (isHailMary ? '#20d8ff' : (isStoneage ? '#b3532e' : (isDisney ? '#1e6fd9' : (isPixar ? '#0b5fff' : (isTelextext ? '#00ffff' : (isFlowerPower1960s ? '#ff4fa1' : (isMarbleGold ? '#c8a24a' : (isRusticWood ? '#8c5a3c' : '#004B44')))))))))))
  var nullColor = isDark ? '#1f2937' : (isHacienda ? '#2f1841' : (isSF1950s ? '#d9cfbf' : (isParis1940s ? '#d4c4ab' : (isHailMary ? '#14243a' : (isStoneage ? '#d9c3a0' : (isDisney ? '#cfe6fb' : (isPixar ? '#d8ecff' : (isTelextext ? '#0c1e2f' : (isFlowerPower1960s ? '#f8e1ff' : (isMarbleGold ? '#e7e6e1' : (isRusticWood ? '#e4ccb0' : '#EFEFEF')))))))))))
  var borderColor = isDark ? 'rgba(148, 163, 184, 0.22)' : (isHacienda ? 'rgba(255, 127, 17, 0.35)' : (isSF1950s ? '#c8b79f' : (isParis1940s ? '#bba78c' : (isHailMary ? 'rgba(125, 249, 255, 0.25)' : (isStoneage ? '#c79c64' : (isDisney ? '#9fc8f8' : (isPixar ? '#8bc5ff' : (isTelextext ? 'rgba(0, 255, 255, 0.25)' : (isFlowerPower1960s ? '#d77bff' : (isMarbleGold ? '#c8b27a' : (isRusticWood ? '#b78f68' : '#e2e8f0')))))))))))
  var titleColor = isDark ? '#ecfdf5' : (isHacienda ? '#fff46b' : (isSF1950s ? '#5a3527' : (isParis1940s ? '#3f2f22' : (isHailMary ? '#c2ff7a' : (isStoneage ? '#4b2f1c' : (isDisney ? '#1e4f9e' : (isPixar ? '#0b5fff' : (isTelextext ? '#ffff00' : (isFlowerPower1960s ? '#9b2fc9' : (isMarbleGold ? '#8a6a23' : (isRusticWood ? '#6f4327' : '#0f172a')))))))))))
  var legendColor = isDark ? '#94a3b8' : (isHacienda ? '#ffb3d7' : (isSF1950s ? '#6f6256' : (isParis1940s ? '#6a5b4c' : (isHailMary ? '#88d6ff' : (isStoneage ? '#6d5845' : (isDisney ? '#3f6ea8' : (isPixar ? '#2f4f78' : (isTelextext ? '#8adfff' : (isFlowerPower1960s ? '#5d3a8f' : (isMarbleGold ? '#6b6560' : (isRusticWood ? '#66503c' : '#64748b')))))))))))
  var axisLineColor = isDark ? '#334155' : (isHacienda ? '#5a2a67' : (isSF1950s ? '#b9a78e' : (isParis1940s ? '#b3a48f' : (isHailMary ? '#2b4f6a' : (isStoneage ? '#b99873' : (isDisney ? '#9fc8f8' : (isPixar ? '#8bc5ff' : (isTelextext ? '#00ffff' : (isFlowerPower1960s ? '#c58cff' : (isMarbleGold ? '#c8b27a' : (isRusticWood ? '#b78f68' : '#cbd5e1')))))))))))
  if (isScifi1950s) {
    plotBg = '#121e34'
    minColor = '#121e34'
    maxColor = '#56f0ff'
    nullColor = '#1d2b45'
    borderColor = 'rgba(86, 240, 255, 0.25)'
    titleColor = '#ffbe5c'
    legendColor = '#8fdfff'
    axisLineColor = '#3eb6d6'
  }
  if (isDrWho) {
    plotBg = '#101b33'
    minColor = '#101b33'
    maxColor = '#4aa8ff'
    nullColor = '#1b2742'
    borderColor = 'rgba(74, 168, 255, 0.25)'
    titleColor = '#7ecbff'
    legendColor = '#a8c8ff'
    axisLineColor = '#4aa8ff'
  }
  if (isPlaySchool) {
    plotBg = '#f2f8ff'
    minColor = '#f2f8ff'
    maxColor = '#4aa8ff'
    nullColor = '#dfeeff'
    borderColor = '#8eb8ff'
    titleColor = '#2f5fb8'
    legendColor = '#4c6492'
    axisLineColor = '#8eb8ff'
  }
  if (isTron) {
    plotBg = '#050a12'
    minColor = '#050a12'
    maxColor = '#00f0ff'
    nullColor = '#0a1828'
    borderColor = 'rgba(0, 240, 255, 0.4)'
    titleColor = '#00f0ff'
    legendColor = '#a8f7ff'
    axisLineColor = '#00f0ff'
  }
  if (isStarWarsANH) {
    plotBg = '#15120b'
    minColor = '#15120b'
    maxColor = '#ffe066'
    nullColor = '#2a2417'
    borderColor = 'rgba(210, 179, 74, 0.28)'
    titleColor = '#ffe066'
    legendColor = '#d9c98b'
    axisLineColor = '#d2b34a'
  }
  if (isAcademic) {
    plotBg = '#f6f9fc'
    minColor = '#f6f9fc'
    maxColor = '#1f4e79'
    nullColor = '#e7eef5'
    borderColor = '#9db4c8'
    titleColor = '#1f4e79'
    legendColor = '#4b6072'
    axisLineColor = '#9db4c8'
  }
  if (isMutedScientist) {
    plotBg = '#f1f5f6'
    minColor = '#f1f5f6'
    maxColor = '#496b7a'
    nullColor = '#e2e9ec'
    borderColor = '#b7c4c9'
    titleColor = '#496b7a'
    legendColor = '#5f6f75'
    axisLineColor = '#b7c4c9'
  }
  if (isTerminal) {
    plotBg = '#0a140a'
    minColor = '#0a140a'
    maxColor = '#33ff66'
    nullColor = '#132613'
    borderColor = 'rgba(51, 255, 102, 0.25)'
    titleColor = '#33ff66'
    legendColor = '#7fe27f'
    axisLineColor = '#33ff66'
  }
  if (isKimiko) {
    plotBg = '#f7f7f4'
    minColor = '#f7f7f4'
    maxColor = '#d9292f'
    nullColor = '#e6e7e5'
    borderColor = '#8ea1b2'
    titleColor = '#0f2a44'
    legendColor = '#4c5f72'
    axisLineColor = '#8ea1b2'
  }
  if (isLondonCalling) {
    plotBg = '#121315'
    minColor = '#121315'
    maxColor = '#e83f7a'
    nullColor = '#1d1f23'
    borderColor = 'rgba(242, 240, 232, 0.22)'
    titleColor = '#e83f7a'
    legendColor = '#cfd4da'
    axisLineColor = '#7f8994'
  }
  if (isGrumpyChemist) {
    plotBg = '#1e2329'
    minColor = '#1e2329'
    maxColor = '#c9a227'
    nullColor = '#2a3038'
    borderColor = 'rgba(201, 162, 39, 0.25)'
    titleColor = '#d4b44a'
    legendColor = '#a8b0ba'
    axisLineColor = '#5c6570'
  }
  if (isOxford) {
    plotBg = '#f4f0e6'
    minColor = '#f4f0e6'
    maxColor = '#002147'
    nullColor = '#e0d8cc'
    borderColor = 'rgba(0, 33, 71, 0.35)'
    titleColor = '#002147'
    legendColor = '#3d3d36'
    axisLineColor = '#002147'
  }
  if (isCambridge) {
    plotBg = '#f0f9ff'
    minColor = '#f0f9ff'
    maxColor = '#0284c7'
    nullColor = '#e0f2fe'
    borderColor = 'rgba(56, 189, 248, 0.35)'
    titleColor = '#0369a1'
    legendColor = '#475569'
    axisLineColor = '#7dd3fc'
  }
  return Highcharts.chart(divName, {
    chart: {
      type: `heatmap`,
      backgroundColor: 'transparent',
      plotBackgroundColor: plotBg
    },
    title: {
      text: chartTitle,
      style: {
        color: titleColor
      }
    },
    colorAxis: {
      min: 0,
      minColor: minColor,
      maxColor: maxColor
    },
    xAxis: {
      title: null,
      lineColor: axisLineColor,
      tickColor: axisLineColor,
      labels: {
        enabled: false
      }

    },
    yAxis: {
      title: null,
      lineColor: axisLineColor,
      tickColor: axisLineColor,
      labels: {
        enabled: false
      }

    },
    legend: {
      layout: `horizontal`,
      itemStyle: {
        color: legendColor
      }
    },
    exporting: {
      enabled: true
    },
    tooltip: (function () {
      let tooltipBg = '#ffffff'
      let tooltipBorder = '#cbd5e1'
      let tooltipText = '#0f172a'
      if (isDark) {
        tooltipBg = '#1e293b'; tooltipBorder = '#334155'; tooltipText = '#f1f5f9'
      } else if (isLondonCalling) {
        tooltipBg = '#121315'; tooltipBorder = '#e83f7a'; tooltipText = '#f2f0e8'
      } else if (isGrumpyChemist) {
        tooltipBg = '#1e2329'; tooltipBorder = '#c9a227'; tooltipText = '#e8e6e0'
      } else if (isOxford) {
        tooltipBg = '#faf6ef'; tooltipBorder = '#002147'; tooltipText = '#1c1b18'
      } else if (isCambridge) {
        tooltipBg = '#ffffff'; tooltipBorder = '#38bdf8'; tooltipText = '#0c4a6e'
      } else if (isKimiko) {
        tooltipBg = '#f7f7f4'; tooltipBorder = '#8ea1b2'; tooltipText = '#1f2b37'
      } else if (isTerminal) {
        tooltipBg = '#0a140a'; tooltipBorder = '#33ff66'; tooltipText = '#c8ffc8'
      } else if (isAcademic) {
        tooltipBg = '#f6f9fc'; tooltipBorder = '#9db4c8'; tooltipText = '#223243'
      } else if (isMutedScientist) {
        tooltipBg = '#f1f5f6'; tooltipBorder = '#b7c4c9'; tooltipText = '#2a3438'
      } else if (isStarWarsANH) {
        tooltipBg = '#15120b'; tooltipBorder = '#d2b34a'; tooltipText = '#f6e8b1'
      } else if (isTron) {
        tooltipBg = '#050a12'; tooltipBorder = '#ff2fd0'; tooltipText = '#e8fdff'
      } else if (isDrWho) {
        tooltipBg = '#101b33'; tooltipBorder = '#4aa8ff'; tooltipText = '#eef7ff'
      } else if (isPlaySchool) {
        tooltipBg = '#f5f9ff'; tooltipBorder = '#8eb8ff'; tooltipText = '#28385a'
      } else if (isScifi1950s) {
        tooltipBg = '#121e34'; tooltipBorder = '#56f0ff'; tooltipText = '#e9fbff'
      } else if (isHacienda) {
        tooltipBg = '#2b1333'; tooltipBorder = '#ff3b8d'; tooltipText = '#fff0f8'
      } else if (isSF1950s) {
        tooltipBg = '#f6efe2'; tooltipBorder = '#c8b79f'; tooltipText = '#2e2a25'
      } else if (isParis1940s) {
        tooltipBg = '#f2e6d2'; tooltipBorder = '#c8b79f'; tooltipText = '#2f2923'
      } else if (isHailMary) {
        tooltipBg = '#0b1322'; tooltipBorder = '#20d8ff'; tooltipText = '#e7fbff'
      } else if (isStoneage) {
        tooltipBg = '#f5e4c4'; tooltipBorder = '#c79c64'; tooltipText = '#2b2218'
      } else if (isDisney) {
        tooltipBg = '#e8f4ff'; tooltipBorder = '#9fc8f8'; tooltipText = '#173055'
      } else if (isPixar) {
        tooltipBg = '#eef7ff'; tooltipBorder = '#8bc5ff'; tooltipText = '#10243f'
      } else if (isTelextext) {
        tooltipBg = '#05101c'; tooltipBorder = '#00ffff'; tooltipText = '#e8faff'
      } else if (isFlowerPower1960s) {
        tooltipBg = '#fff2ff'; tooltipBorder = '#d77bff'; tooltipText = '#3b235a'
      } else if (isMarbleGold) {
        tooltipBg = '#f7f6f2'; tooltipBorder = '#c8b27a'; tooltipText = '#2f3136'
      } else if (isRusticWood) {
        tooltipBg = '#f5e7d2'; tooltipBorder = '#b78f68'; tooltipText = '#2e2218'
      }
      return {
        backgroundColor: tooltipBg,
        borderColor: tooltipBorder,
        style: { color: tooltipText }
      }
    })(),
    series: {
      nullColor: nullColor,
      borderColor: borderColor,
      borderWidth: 1,
      type: `heatmap`
    }

  })
}

/**
 * Create a HighCharts pie chart. These are used on the basecalled data tab.
 * @param divId {string} The id of the div to render chart to.
 * @param chartTitle {string} The title of the chart.
 * @return {*}
 */
function makePieChart (divId, chartTitle) {
  return Highcharts.chart(divId, {
    chart: {
      plotBackgroundColor: null,
      plotBorderWidth: null,
      plotShadow: false,
      type: `pie`
    },
    title: {
      text: chartTitle
    },
    tooltip: {
      pointFormat: `{series.name}: <b>{point.percentage:.1f}%</b>`
    },
    accessibility: {
      point: {
        valueSuffix: `%`
      }
    },
    lang: {
      noData: `Looking for data probably.`
    },
    plotOptions: {
      pie: {
        allowPointSelect: true,
        cursor: `pointer`,
        dataLabels: {
          enabled: true,
          format: `<b>{point.name}</b>: {point.percentage:.1f} %`
        }
      }
    },
    series: []
  })
}

/**
 * Add the tooltips to help icons
 */
function addJqueryTooltip () {
  $(`.masterTooltip`).on(`mouseover`, (event) => {
    const text = event.target.attributes.tooltip.textContent
    $(`body`)
      .append(`<div class="toolTip" style="left:${event.pageX + 25}px;
 top: ${event.pageY - 25}px; display: inline-block">${text}</div>`)
  }).on(`mouseout`, event => {
    $(`.toolTip`).remove()
  })
}

/**
 * Transform yield in bases to make human readable yields for table columns i.e 100Mb, 1Gb etc.
 * @param data {[]} Data for the whole table
 * @param type {string} Data for the
 * @param row  {[]} data for the row
 * @return {string|*}
 */
function humanReadableYield (data, type, row) {
  // The below function returns the yield in a human readable format
  if (type === `display` || type === `filter`) {
    var UNITS = [``, `k`, `M`, `G`, `T`, `P`, `E`, `Z`]
    var factor = 1000
    var suffix = `b`
    for (var i = 0; i < UNITS.length; i++) {
      if (Math.abs(data) < factor) {
        data = data.toFixed(2)
        if (/\.00$/.test(data)) {
          data = data.substr(0, data.length - 3)
        }
        return `${`${data} ${UNITS[i]}${suffix}`}`
      }
      data /= factor
    }
    data = data.toFixed(2).replace(/\B(?=(\d{3})+(?!\d))/g, `,`)
    if (/\.00$/.test(data)) {
      data = data.substr(0, data.length - 3)
    }
    return `<a href="/web/private/${that._linkDestination}/${row.id}/">${`${data} Y${suffix}`}</a>`
  } else {
    return data
  }
}
