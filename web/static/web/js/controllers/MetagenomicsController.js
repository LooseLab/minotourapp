class MetagenomicsController {
  /**
   *
   * @param flowcellId {number} The primary key of the flowcell record in the database
   */
  constructor (flowcellId) {
    /**
     */
    this._flowcellId = flowcellId
    this._axiosInstance = axios.create({
      headers: { "X-CSRFToken": getCookie(`csrftoken`) }
    })
    setSelectedBarcode(`All reads`, `Metagenomics`)
    this._selectedBarcode = getSelectedBarcode(`Metagenomics`)
    this._addTooltips()
    this._addExpandButtonListener()
    this.colour = d3.scaleOrdinal(d3.schemeCategory10)
    this._first = true
    this._range = $(`.input-range`)
    this._value = $(`.taxa-level`)
    this._taxas = [`species`, `genus`, `family`, `order`, `class`, `phylum`, `superkingdom`]
    // the taxa titles we wish to display under the slider
    this._displayTaxas = [`Species`, `Genus`, `Family`, `Order`, `Class`, `Phylum`, `Kingdom`]
    this._donutData = []
    this._barcodesSet = new Set()
    this._barcodeDivElement = $(`#nav-tabs-meta-barcodes`)
    this._addListenerToRangeSlider()
  }

  updateTab () {
    this._getTotalReadsTable(this._flowcellId)
    this._drawDonutAndTable(this._flowcellId)
    this._metaDataTable(this._flowcellId)
    this._addBarcodeTabs(this._flowcellId)
    this._updateMappingTable(this._flowcellId)
    this._drawSimpleTable(this._flowcellId)
    this._initialiseSankey(this._flowcellId)
  }

  _addExpandButtonListener () {
    // Changes the hide data button text to show/hide based on click
    $(`#expand-button`).click(function () {
      $(this).text((i, old) => old.replace(/(\r\n|\n|\r|\s)/gm, ``) === `Expanddata` ? `Hide data` : `Expand data`)
    })
  }

  _addTooltips () {
    // function actually totally unrelated to the donut, does the tooltips for the help icons
    // d3.selectAll(`.masterTooltip`).on(`mouseover`, function (event, d) {
    //   const text = d3.select(this).attr(`tooltip`)
    //   const tooltip = d3.select(`body`)
    //     .append(`div`)
    //     .attr(`class`, `toolTip`)
    //   tooltip
    //     .style(`left`, `${event.pageX + 25}px`)
    //     .style(`top`, `${event.pageY - 25}px`)
    //     .style(`display`, `inline-block`)
    //     .html(text)
    // }).on(`mouseout`, d => {
    //   d3.select(`.toolTip`).remove(`*`)
    // })
  }

  _revealMetagenomicsPage () {
    d3.select(`#loading-sign`).transition().duration(1000).style(`opacity`, 0)
    setTimeout(() => {
      $(`body`).addClass(`loaded`)
      d3.select(`#loading-sign`).style(`display`, `none`)
      d3.select(`.vis-container`).style(`display`, `contents`)
    }, 1500)
  }

  _getTotalReadsTable (flowcellId) {
    // Get total reads table updates the total reads table at te bottom of the page
    // Get data from the api
    // Jquery selector
    const that = this
    const table = $(`.tableLand`)
    // If the table already exists, use the DataTable APi to update in place
    if ($.fn.DataTable.isDataTable(table)) {
      table.DataTable().ajax.reload(null, false)
    } else {
      // else the databale must be initialised
      table.DataTable({
        initComplete: (settings, json) => {
          that._revealMetagenomicsPage()
        },
        processing: true,
        serverSide: true,
        ajax: {
          url: `/api/v1/metagenomics/${flowcellId}/table/`,
          async: true,
          error: (xhr, error, code) => {
            console.error(xhr)
            console.error(code)
          }
        },
        columns: [
          { data: `barcode_name` },
          { data: `superkingdom` },
          { data: `phylum` },
          { data: `classy` },
          { data: `order` },
          { data: `family` },
          { data: `genus` },
          { data: `species` },
          { data: `num_matches` },
          { data: `proportion_of_classified` }
        ]
      }
      )
    }
  }

  move (event, d) {
    console.log(event)
    console.log(d)
    d3.select(`.${event.sourceEvent.srcElement.firstChild.classList[0]}`).attr(`transform`, event.transform)
  }

  _drawPie (countedData, pie, arc, svg) {
    // draw the donut chart
    // select an actual good colour scheme
    const color = d3.scaleOrdinal(d3.schemeCategory10)
    // Select all the slices and bind the provided data after it's been transformed by d3.pie()
    const slice = svg.select(`.slices`).selectAll(`path.slice`)
      .data(pie(countedData))
    // enter your data, returned from pie(countedData, insert a path, set d to data provided by arc function)
    // Update existing donut slices
    slice.attr(`class`, `slice`)
      .style(`fill`, (d, i) => color(i))
      .attr(`stroke`, `black`)
      .attr(`stroke-width`, 0.2)
      .attr(`d`, arc).select(`title`)
      .text(d => `${d.data.label}\n${d.value} reads`)
    // add the new donut slices, held in d3s enter selection
    slice.enter()
      .insert(`path`)
      .attr(`class`, `slice`)
      .style(`fill`, (d, i) => color(d.data.label))
      .attr(`stroke`, `black`)
      .attr(`stroke-width`, 0.2)
      .attr(`d`, arc).append(`title`)
      .text(d => `${d.data.label}\n${d.value} reads`)
    // remove any slices held in the exit selection, that used to have DOM elements but now have no data for them
    slice.exit().remove()
  }

  _drawDonutAndTable (flowcellId) {
    const selectedBarcode = getSelectedBarcode(`Metagenomics`)
    const rankTable = d3.select(`.rank-table-container`)
    // setup the donut chart
    // the taxas in the order we want, to access them from the AJAX get request results
    // Calculate the width
    const width = ($(window).width() * 0.25) - 76
    // Calculate the height
    const height = $(window).height() * 0.27
    // the radius, the smallest of the width and height /2 so it fits in hte svg
    const radius = Math.min(width, height) / 2
    // The d3 zoom function
    const zoom = d3.zoom()
      .scaleExtent([1, 10]).translateExtent([[0, 0], [width, height]])
      .on(`zoom`, this.move)
    // Select the Range of the html slider and the Displayed label underneath the slider
    const value = $(`.taxa-level`)
    // Declare the svg and group element variables
    let svg
    let g
    // pie is a d3 function that transforms the value for each of the .values in the array of objects, into
    const pie = d3.pie().padAngle(0.01).sort(null)
      .value(d => d.value)
    // arc is a d3 function that generates arc paths from the data provided from d3.pie
    const arc = d3.arc()
    // \The inner and outer radius of the actual donut slices
      .innerRadius(radius * 0.8)
      .outerRadius(radius * 0.5)
      .padAngle(0.02)
    this._pie = pie
    this._arc = arc
    // If there is already a donut-svg
    if ($(`.donut-svg`).length !== 0) {
      // Select the svg
      svg = d3.select(`.donut-svg`)
      // recenter, as resizing doesn't change the transformation
      d3.select(`.slices`).attr(`transform`, `translate(${width / 2},${height / 2})`)
    } else {
      // If no extant svg, append a new svg and g element
      svg = d3.select(`.donutContainer`).append(`svg`)
        .attr(`class`, `donut-svg`)
        .attr(`width`, width)
        .attr(`height`, height)
        .attr(`margin`, `auto`)
        .attr(`display`, `block`)
        .call(zoom)
      this._svg = svg
      g = svg.append(`g`).attr(`class`, `badCopNoDonut`)
      // append g elements to the svg, attr transform translate centers them in the svg, otherwise drawn offscreen
      g.append(`g`)
        .attr(`class`, `slices`)
        .attr(`transform`, `translate(${width / 2},${height / 2})`)
      // Set the slider level to Display Species
      value.html(this._displayTaxas[0])
    }
    // Get the data from the server to sisplay
    this._axiosInstance.get(`/api/v1/metagenomics/donut`, {
      params: {
        flowcellId,
        visType: `donut`,
        barcode: selectedBarcode
      }
    }).then(
      result => {
        // if there is no data return and try again when interval is up on $interval
        if (result.data === undefined || result.data.species === undefined) {
          svg.select(`.badCopNoDonut`).selectAll(`*`).remove()
          svg.append(`text`).text(`No data to display`).attr(`text-anchor`, `middle`).attr(`x`, `50%`).attr(`y`, `50%`)
          // Now the same for the table
          rankTable.selectAll(`*`).remove()
          rankTable.append(`div`).style(`height`, `${height}px`).style(`transform`, `translateY( 50%)`).attr(`class`, `col-md-12`).append(`text`).text(`No data to display`).classed(`no-data`, true)
          return
        } else {
          rankTable.select(`div`).remove()
        }
        this._donutData = result.data
        const dataToDraw = result.data
        // get the right taxa string, so we can use it as a key o the results object
        const currentlySelectedTaxa = this._taxas[this._range.val()]
        // data1 is the data for the currently selected taxa from the results array
        const data1 = dataToDraw[currentlySelectedTaxa]
        // Draw a new pie chart and table
        this._drawPie(data1, pie, arc, svg)
        const tableData = this._formatTableData(data1, this.colour)
        this._drawTables(`rank-table-container`, tableData)
      }).catch(
      error => {
        console.error(error)
      }
    )
  }

  _addListenerToRangeSlider () {
    /**
     * Add the on change listener to redraw the charts to the range slider.
     */
    const that = this
    this._range.on(`input`, () => {
      // the selected number for the slider level (0-6)
      const number = that._range.val()
      // get the right taxa for the key to the results object
      const currentSelectedTaxa = that._taxas[number]
      // get the results data array for the currently selected taxa clade
      const sortedData = this._donutData[currentSelectedTaxa]
      // set the html below the slider to the right level
      that._value.html(that._displayTaxas[number])
      // datalength - how many members ar ein this clade (1-20)
      // draw a new donut
      that._drawPie(sortedData, that._pie, that._arc, that._svg)
      const tableData = that._formatTableData(sortedData, that.colour)
      this._drawTables(`rank-table-container`, tableData)
    })
  }

  _formatTableData (countedData, color) {
    /**
     * @type {*[]}
     */
    const tableData = []
    countedData.forEach((datum, index) => {
      tableData.push({
        Species: datum.label,
        "# Matches": datum.value,
        Rank: index + 1,
        Key: color(index)
      })
    })
    return tableData
  }

  _drawTables (selection, dataToDraw) {
    // draw the table using d3
    // get table height
    const columns = [`Species`, `# Matches`, `Rank`, `Key`]
    let rows
    let cells
    const table = this._first ? d3.select(`.${selection}`).append(`table`).attr(`class`, `table table-hover taxons`) : d3.select(`.${selection}`).select(`table`)
    const tbody = this._first ? table.append(`tbody`).attr(`class`, `${selection}Body`) : table.select(`tbody`)
    if (this._first) {
      table.append(`thead`).append(`tr`)
        .selectAll(`th`)
        .data(columns).enter()
        .append(`th`)
        .text(column => column)
    }
    rows = tbody.selectAll(`tr`)
    rows.data(dataToDraw).exit().remove()
    rows.data(dataToDraw).enter().append(`tr`)
    rows.data(dataToDraw).exit().remove()
    rows = tbody.selectAll(`tr`)
    cells = rows.selectAll(`td`)
    cells.data(row => columns.map(column => {
      // if (column !== "Key"){
      return { column: column, value: row[column] }
      // }
    })).enter().append(`td`)
    cells = rows.selectAll(`td`)
    cells.data(row => columns.map(column => ({ column, value: row[column] }))).attr(`class`, (d, i) => columns[i])
      .style(`background-color`, (d, i) => {
      // if the cell contains the key, set the background colour to the rgb value in the data
        if (typeof d.value === `string`) {
          if (d.value.includes(`#`)) {
            return d.value
          } else {
          // else set the background to white
            return `white`
          }
        }
      })
    // fill the cell with the string by setting the inner html
      .html(d => {
        if (typeof d.value === `string`) {
          if (!d.value.includes(`#`)) {
            return d.value
          }
        } else {
          return d.value
        }
      })
    this._first = false
  }

  _metaDataTable (flowcellId) {
    // AJAX request for the metadata
    this._axiosInstance.get(`/api/v1/metagenomics/metagenomics-metadata/`, { params: { flowcellId } }).then(
      result => {
        // Metaheader draws or updates the metadata header at the top of the visualisation page
        // initialise variables
        let table, head, row
        const data = result.data.result
        // if we have validation results for this flowcell
        const targetDiv = result.data.validation ? `.metadata-table` : `.metadata-table-alt`
        if (!d3.select(targetDiv).classed(`has-table`)) {
          // append the table to the page
          table = d3.select(targetDiv).classed(`has-table`, true).append(`table`).attr(`class`, `table table-striped`).attr(`table-layout`, `fixed`)
          head = table.append(`thead`)
          row = head.append(`tr`)
        } else {
          // select the already existing row
          row = d3.select(targetDiv).select(`table`).select(`tr`)
        }
        // If there is no data, return
        if (data[1].value === 0) {
          return
        }
        // Add any DOM elements for data in the enter selection, which is if there is no th cells already
        row.selectAll(`th`).data(data).enter().append(`th`)
        // Update the innerHTML value of the cells
        row.selectAll(`th`).data(data).html(d => d.key + d.value)
      }
    ).catch(
      error => {
        console.error(
          error
        )
      })
  }

  _addBarcodeTabs (flowcellId) {
    this._axiosInstance.get(`/api/v1/metagenomics/${flowcellId}/metagenomic-barcodes`).then(
      response => {
        const barcodes = response.data.data.sort()
        const tabs = response.data.tabs
        const selectedBarcode = getSelectedBarcode(`Metagenomics`)
        const barcodeTabElements = []
        const alertLevels = { 0: `green-alert-tab`, 1: `green-alert-tab`, 2: `orange-alert-tab`, 3: `red-alert-tab` }
        // add new fetched barcode tabs to the set
        barcodes.forEach(barcode => {
          this._barcodesSet.add(barcode)
        })
        this._barcodesSet.forEach(barcode => {
          let tabLevel = ``
          let active = ``
          if (tabs.length) {
            tabLevel = alertLevels[tabs[barcode]]
          }
          if (barcode === selectedBarcode) {
            active = `active`
          }
          barcodeTabElements.push(`<li class="barcode-meta-tab nav-item ${tabLevel}"><a class="nav-link ${active}" id="meta-barcode-link">${barcode}</a></li>`)
        })
        this._barcodeDivElement.html(barcodeTabElements)
        $(`#meta-barcode-link`).on(`click`, () => {
          const clickedBarcode = event.srcElement.innerText
          setSelectedBarcode(clickedBarcode, `Metagenomics`)
          this.updateTab()
        })
      }
    ).catch(
      error => {
        console.error(error)
      }
    )
  }

  _compare (a, b) {
    if (a.Species < b.Species) {
      return -1
    }
    if (a.Species > b.Species) {
      return 1
    }
    return 0
  }

  _compare2 (a, b) {
    // sorting function to put results in alphabetical order
    if (a[`Validation species`] < b[`Validation species`]) { return -1 }
    if (a[`Validation species`] > b[`Validation species`]) { return 1 }
    return 0
  }

  _updateMappingTable (flowcellId) {
    let table
    let thead
    let tbody
    let rows
    let cells
    let rowCount = 0
    const columns = [`Alert level`, `Species`, `Num. matches`, `Prop. classified (%)`,
      `Num. mapped`,
      `Mapped prop. matches (%)`, `Target reads`,
      `Target prop. mapped (%)`
    ]
    const barcode = getSelectedBarcode(`Metagenomics`)
    // Order the results correctly
    this._axiosInstance.get(`/api/v1/metagenomics/mapped-targets`, { params: { flowcellId, barcode } }).then(
      result => {
        const tableData = result.data.table
        const tableInitialised = d3.select(`.alert-table`).classed(`has-tabley?`)
        if (result.status === 204) {
          d3.select(`.alert-table-complex`).remove()
          d3.select(`.button-row`).remove()
          d3.select(`#demo`).classed(`in`, true)
          return
        }
        if (tableData === undefined) {
          const alertTable = d3.select(`.alert-table`)
          alertTable.selectAll(`*`).remove()
          alertTable.append(`text`).text(`No data to display`)
          alertTable.classed(`no-data`, true)
          return
        } else {
          d3.select(`.alert-table`).classed(`no-data`, false).select(`text`).remove()
        }
        tableData.sort(this._compare)
        table = tableInitialised ? d3.select(`.alert-table`).select(`table`) : d3.select(`.alert-table`).classed(`has-tabley?`, true).style(`width`, `100%`).append(`table`).attr(`class`, `table map-alert`)
        thead = tableInitialised ? table.select(`thead`) : table.append(`thead`).append(`tr`)
        tbody = tableInitialised ? table.select(`tbody`) : table.append(`tbody`).attr(`class`, `alert-tbody`)
        thead.selectAll(`th`)
          .data(columns).enter()
          .append(`th`)
          .text(column => column)
        rows = tbody.selectAll(`tr`)
        rows
          .data(tableData)
          .enter()
          .append(`tr`).attr(`id`, d => d.Species.replace(/ /g, `_`))
        rows.data(tableData).exit().remove()
        rows = tbody.selectAll(`tr`)
        cells = rows.selectAll(`td`)
        cells.data(row => columns.map(column => {
          return { column, value: row[column] }
        }))
          .enter()
          .append(`td`)
        cells = rows.selectAll(`td`)
        cells.data(row => columns.map(column => {
          return { column: column, value: row[column] }
        })).attr(`id`, (d, i) => {
          if (d.column === `Alert level`) {
            rowCount += 1
          }
          return `${columns[i].replace(` `, `_`)}_${rowCount.toString()}`
        })
          .style(`background-color`, function (d, i) {
          // if the cell contains the key, set the background colour to the rgb value in the data
            const variable = d3.select(`#${d3.select(this).node().parentNode.children[0].id}`)
            if (d.column === `Num. mapped` && d.value > 0) {
              variable.classed(`green-alert`, false)
              variable.classed(`yellow-alert`, false)
              variable.classed(`red-alert`, true)
            } else if (d.column === `Target reads` && d.value > 0) {
              variable.classed(`green-alert`, false)
              variable.classed(`orange-alert`, false)
              variable.classed(`yellow-alert`, false)
              variable.classed(`red-alert`, true)
            } else if (d.column === `Num. matches` && d.value > 0) {
              variable.attr(`class`, `green-alert`)
            } else if (d.column === `Num. matches` && d.value === 0) {
              variable.attr(`class`, `green-alert`)
            }
          })
        // fill the cell with the string by setting the inner html
          .html(d => d.value)
      })
  }

  _drawSimpleTable (flowcellId) {
  // declare variables at scope start
    let table; let thead; let tbody; let rows; let cells; let data; let limit
    const columns = [`Status`, `Low probability`, `Validation species`, `Detected`]
    const barcode = getSelectedBarcode(`Metagenomics`)
    // get the alerts
    this._axiosInstance.get(`/api/v1/metagenomics/alerts`, { params: { flowcellId, barcode } }).then(
      response => {
        // If there is no validation mapping results objects in the database for this run
        if (response.status === 204) {
          d3.select(`#validation`).remove()
          d3.select(`#analysis`).classed(`show`, true)
          return
        }
        data = response.data.table
        limit = response.data.conf_detect_limit
        // sort the data aphabetically on the Validation species column
        data.sort(this._compare2)
        // if the table has been created on the page drawing, select it so we can update the existing table
        const tableInitialised = d3.select(`.simple-alert-table`).classed(`has-tabley?`)
        table = tableInitialised ? d3.select(`.simple-alert-table`).select(`table`) : d3.select(`.simple-alert-table`).classed(`has-tabley?`, true).style(`width`, `100%`).append(`table`).attr(`class`, `table map-alert`)
        thead = tableInitialised ? table.select(`thead`) : table.append(`thead`).append(`tr`)
        tbody = tableInitialised ? table.select(`tbody`) : table.append(`tbody`).attr(`class`, `alert-tbody`)
        // select all the th elements
        thead.selectAll(`th`)
        // bind the columns array
          .data(columns).enter()
        // enter the selection, so enter is a selection of data that isn't already bound to an html element
          .append(`th`).style(`text-align`, `center`)
        // give it an id matching the column name
          .attr(`id`, column => column.replace(` `, `_`)).attr(`class`, (d, i) => {
            // give it a class of it's index base 0 - used for to populate the table below
            return i
            // make the status column 10% of the table width
          }).style(`width`, d => {
            if (d === `Status`) {
              return `10%`
            }
            // return the column text
          }).text(column => column)

        // select all the rows
        rows = tbody.selectAll(`tr`)
        // same binding of new elements as above
        rows
          .data(data)
          .enter()
          .append(`tr`).attr(`id`, d => {
            // give each row an id of the species title
            return d[`Validation species`].replace(/ /g, `_`)
          })
        // exit is a selection of html elements who's bound data in the array has been removed, so remove the html
        rows.data(data).exit().remove()
        // select all the new rows
        rows = tbody.selectAll(`tr`)
        // select all the cells in the rows
        cells = rows.selectAll(`td`)
        // for each column in each row, pass a dictionary down the chain of column and that columns value in each row
        cells.data(row => columns.map(column => ({ column: column, value: row[column] })))
        // enter the new data, anything not bound to a html element append a cell
          .enter()
          .append(`td`)
        // select all the newly created and exisiting html td elements
        cells = rows.selectAll(`td`)
        // for each row return a mapped dictionary of all possible values for that rows columns
        cells.data(row => {
          return columns.map(column => ({
            column,
            value: row[column],
            read_count: row.read_count,
            detected: row.Detected,
            species: row[`Validation species`],
            conf_limit: row.conf_limit,
            detected_at: row.detected_at
          }))
          // set the id of the cell as the column id
        }).attr(`id`, (d, i) => columns[i].replace(` `, `_`)).style(`background-color`, d => {
          // change status cell colour depending on findings
          if (d.column === `Status` && d.detected === false) {
            return `green`
          } else if (d.column === `Status` && d.detected === true) {
            // a nice pastely red colour
            return `rgba(255, 0, 0, 0.7)`
          }
        }).html(
          d => {
            // set the validation species html to not have a underscore
            d3.select(`#Validation species`).html(`Validation species`)
            // the code below places the html in the left cell if we have seen lots of reads but no matches
            if (d.column === `Low probability` && d.detected === false && d.conf_limit >= 50000) {
              d3.select(`#Low_probability`).html(`Low risk, less than 1 in ${limit} reads`)
              return d.species
            } else if (d.column === `Validation species` && d.detected === false && d.conf_limit < 50000) {
              d3.select(`#Potential_threats`).html(`Validation species (< 1 in ${limit} reads)`)
              return d.species
            } else if (d.column === `Detected` && d.detected === true) {
              return `${d.species} (1 in ${d.detected_at} reads)`
            }
          }
        )
      })
  }

  _drawSankeyDiagram (nodesObj, sankey, g, format, color, width) {
  // Draw the diagram
  // use d3-sankey to sankeyify the data, providing path coordinates
    const link = g.append(`g`)
    sankey(nodesObj)
    // draw the links
    link.attr(`class`, `links`)
      .attr(`fill`, `none`)
      .attr(`stroke-opacity`, 0.5)
      .selectAll(`path`)
      .data(nodesObj.links)
  .join("path")
    .attr("d", d3.sankeyLinkHorizontal())
    .attr(`stroke`, (d, i) => color(d.path))
      .attr(`stroke-width`, d => Math.max(0.5, d.width)).append(`title`)
      .text(d => `${d.source.name} → ${d.target.name}\n${format(d.value)}`)
    // append the title
    // link.append(`title`)
    //   .text(d => `${d.source.name} → ${d.target.name}\n${format(d.value)}`)
    // Append the nodes to the svg
    g.append(`g`).attr(`class`, `nodes`)
      .selectAll(`rect`)
      .data(nodesObj.nodes)
      .enter()
      .append(`rect`)
      .attr(`x`, d => d.x0)
      .attr(`y`, d => d.y0)
      .attr(`height`, d => d.y1 - d.y0)
      .attr(`width`, d => d.x1 - d.x0)
      .attr(`fill`, d => color(d.path))
      .append(`title`)
      .text(d => `${d.name}\n${format(d.value)}`)
    // append the text labels to the svg
    g.append(`g`).attr(`class`, `text`)
      .style(`font`, `10px sans-serif`)
      .selectAll(`text`)
      .data(nodesObj.nodes)
      .enter().append(`text`)
      .attr(`x`, d => d.x0 < width / 2 ? d.x1 + 6 : d.x0 - 6)
      .attr(`y`, d => (d.y1 + d.y0) / 2)
      .attr(`dy`, `0.35em`)
      .attr(`text-anchor`, d => d.x0 < width / 2 ? `start` : `end`)
      .text(d => d.name)
  }

  // update or draw the existing svg using the AJAX results from the server
  _updateSankey (flowcellId, sankey, checkForData, svg, g, format, color, width, selectedBarcode) {
  // TODO species limit one day
    this._axiosInstance.get(`/api/v1/metagenomics/sankey`, { params: { flowcellId, barcode: selectedBarcode } }).then(
      result => {
        const sankeyData = result.data.sankey; const hasRun = result.data.run
        console.log(hasRun)
        // if theres no data from the server
        if (sankeyData === undefined && hasRun === true) {
          svg.select(`.contain`).selectAll(`*`).remove()
          svg.append(`text`).text(`No data to display`).attr(`text-anchor`, `middle`).attr(`x`, `50%`).attr(`y`, `50%`)
          return
        } else if (hasRun === false) {
          return
        } else if (hasRun === true) {
          d3.select(`.sankeyContainer`).style(`display`, `block`)
        }
        // If there is data and checkForData is true, clear the css elements and loading sign so we can see the graphics
        // TODO update in place
        svg.select(`.contain`).selectAll(`*`).remove()
        this._drawSankeyDiagram(sankeyData, sankey, g, format, color, width)
      }).catch(
      error => {
        console.error(error)
      })
  }

  // top level function
  _initialiseSankey (flowcellId) {
    const selectedBarcode = getSelectedBarcode(`Metagenomics`)
    // Check the tab value
    // set svg width and height
    // height of page
    let hi = $(window).height() * 0.6
    // width of page
    const width = $(window).width() - 120
    let svg
    let g
    const formatNumber = d3.format(`,.0f`)
    const format = d => `${formatNumber(d)} Reads`
    // create colour scheme
    const color = d3.scaleOrdinal(d3.schemeCategory10)
    const checkForData = true
    hi = d3.max([hi, 480])
    const sankey = d3.sankey()
      .nodeWidth(20)
      .nodePadding(3)
      .size([width, hi * 0.95]).nodeId(d => d.name)
    // panning and zoom function

    if ($(`.svg-sankey`).length) {
      svg = d3.select(`.svg-sankey`)
      g = d3.select(`.contain`)
    } else {
      svg = d3.select(`.svg-container`).append(`svg`).attr(`width`, width).attr(`class`, `svg-sankey`)
        .attr(`height`, hi)
      g = svg.append(`g`).attr(`class`, `contain`)
    }
    this._updateSankey(flowcellId, sankey, checkForData, svg, g, format, color, width, selectedBarcode)
  }
}
