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
    this._addListenerToRangeSlider()
  }

  updateTab () {
    this._getTotalReadsTable(this._flowcellId)
    this._drawDonutAndTable(this._flowcellId)
    this._metaDataTable(this._flowcellId)
  }

  _addExpandButtonListener () {
    // Changes the hide data button text to show/hide based on click
    $(`#expand-button`).click(function () {
      $(this).text((i, old) => old.replace(/(\r\n|\n|\r|\s)/gm, ``) === `Expanddata` ? `Hide data` : `Expand data`)
    })
  }

  _addTooltips () {
    // function actually totally unrelated to the donut, does the tooltips for the help icons
    d3.selectAll(`.masterTooltip`).on(`mouseover`, function () {
      const text = d3.select(this).attr(`tooltip`)
      const tooltip = d3.select(`body`)
        .append(`div`)
        .attr(`class`, `toolTip`)
      tooltip
        .style(`left`, d3.event.pageX + 25 + `px`)
        .style(`top`, d3.event.pageY - 25 + `px`)
        .style(`display`, `inline-block`)
        .html(text)
    }).on(`mouseout`, function (d) {
      d3.select(`.toolTip`).remove(`*`)
    })
  }

  _revealMetagenomicsPage () {
    d3.select(`#loading-sign`).transition().duration(1000).style(`opacity`, 0)
    setTimeout(function () {
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

  move () {
    d3.select(`.badCopNoDonut`).attr(`transform`, d3.event.transform)
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
    const width = ($(window).width() * 0.25) - 50
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
      .value(function (d) {
        return d.value
      })
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
      d3.select(`.slices`).attr(`transform`, `translate(` + width / 2 + `,` + height / 2 + `)`)
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
        .attr(`transform`, `translate(` + width / 2 + `,` + height / 2 + `)`)
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
    const table = this._first ? d3.select(`.${selection}`).append(`table`).attr(`class`, `table table-hover taxons`) : d3.select(`.` + selection).select(`table`)
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
    cells.data(row => columns.map(column => ({ column: column, value: row[column] }))).attr(`class`, function (d, i) {
      return columns[i]
    })
      .style(`background-color`, (d, i) => {
      // if the cell contains the key, set the background colour to the rgb value in the data
        if (typeof d.value === `string`) {
          if (d.value.includes(`#`)) { return d.value } else {
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
    this._axiosInstance.get(`/api/v1/metagenomics/centrifuge_metadata/`, { params: { flowcellId } }).then(
      result => {
        // Metaheader draws or updates the metadata header at the top of the visualisation page
        // initialise variables
        let table, head, row
        const data = result.data.result
        const targetDiv = result.data.validation ? `.meta_taberu` : `.meta_taberu_alt`
        if (this._first) {
          // append the table to the page
          table = d3.select(targetDiv).append(`table`).attr(`class`, `table table-striped`).attr(`table-layout`, `fixed`)
          head = table.append(`thead`)
          row = head.append(`tr`)
        } else {
          // select the already existing row
          row = d3.select(targetDiv).select(`table`).select(`tr`)
        }
        // If there is no data, return
        if (data[1].value === 0) { return }
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
}
