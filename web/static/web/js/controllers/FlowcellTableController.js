/**
 * Controller for the tables on the Flowcell page and the Flowcell Manager page.
 */
class FlowcellTableController {
    /**
     *
     * @param linkDestination {string} The destination of the formatted link. flowcell_manager or flowcells.
     * @param elementId {string} The target Element to render this table inside.
     */
    constructor(linkDestination, elementId) {
        this._elementId = elementId;
        this._linkDestination = linkDestination;
        this._datatableObj = null;
        this._interval = setInterval(this.refreshTable, 30000);
        // ON page unload clear the interval, so when user navigates away from table
        $(window).on("unload", function () {
            console.log("stopping refresh");
            clearInterval(this._interval)
        })
    }

    /**
     * @function refresh the flowcell table.
     */
    refreshTable() {
        console.log("Starting refresh");
        this._datatableObj.DataTable().ajax.reload(null, false)
    }

    /**
     * @function Draw the flowcells index table.
     */
    renderTable() {
        let that = this;
        var table_all_runs;

        $.fn.dataTable.render.format_link = function () {
            return function (data, type, full, meta) {

                return `<a href="/web/private/${that._linkDestination}/${full["id"]}/">${data}</a>`;
            };
        };

        this._datatableObj = $(`#${this._elementId}`).DataTable({
            // callback for each created row, add styling class so cursor is pointer showing it's a link
            "createdRow": function (row, data, index) {
                $(row).addClass("stylable");
            },

            "scrollX": true,
            "order": [[3, "desc"]],
            'ajax': {
                'url': "/api/v1/flowcells/",

            },

            'columnDefs': [
                {
                    'targets': 0,
                    'data': "name",
                    'render': $.fn.dataTable.render.format_link()

                },
                {
                    'targets': 1,
                    'data': "size",
                    'render': $.fn.dataTable.render.format_link()
                },
                {
                    'targets': 2,
                    'data': "sample_name",
                    'render': $.fn.dataTable.render.format_link()
                },
                {
                    'targets': 3,
                    'data': "start_time",
                    'render': function (data, type, full, meta) {
                        return `<a href="/web/private/${that._linkDestination}/${full["id"]}/">${moment(data).format('YYYY-MM-DD HH:mm')}</a>`;
                    }
                },
                {
                    'targets': 4,
                    'data': "number_reads",
                    'render': $.fn.dataTable.render.format_link()
                },
                {
                    'targets': 5,
                    'data': "number_reads_processed",
                    'render': $.fn.dataTable.render.format_link()
                },
                {
                    'targets': 6,
                    'data': "number_runs",
                    'render': $.fn.dataTable.render.format_link()
                },
                {
                    'targets': 7,
                    'data': "number_barcodes",
                    'render': $.fn.dataTable.render.format_link()
                },
                {
                    'targets': 8,
                    'data': "total_read_length",
                    render: function (data, type, row) {
                        // The below function returns the read length in a human readable format
                        if (type === "display" || type === "filter") {
                            var UNITS = ["", "k", "M", "G", "T", "P", "E", "Z"];
                            var factor = 1000;
                            var suffix = "b";
                            for (var i = 0; i < UNITS.length; i++) {
                                if (Math.abs(data) < factor) {
                                    data = data.toFixed(2);
                                    if (/\.00$/.test(data)) {
                                        data = data.substr(0, data.length - 3);
                                    }
                                    return `<a href="/web/private/${that._linkDestination}/${row["id"]}/">${`${data} ${UNITS[i]}${suffix}`}</a>`;
                                }
                                data /= factor;
                            }
                            data = data.toFixed(2).replace(/\B(?=(\d{3})+(?!\d))/g, ",");
                            if (/\.00$/.test(data))
                                data = data.substr(0, data.length - 3);
                            return `<a href="/web/private/${that._linkDestination}/${row["id"]}/">${`${data} Y${suffix}`}</a>`;
                        } else {
                            return data
                        }
                    }
                },
                {
                    'targets': 9,
                    'data': "average_read_length",
                    'render': $.fn.dataTable.render.format_link()
                },
                {
                    'targets': 10,
                    'data': "is_active",
                    'render': $.fn.dataTable.render.format_link()
                },
                {
                    'targets': 11,
                    'data': "has_fastq",
                    'render': $.fn.dataTable.render.format_link()
                },
                {
                    'targets': 12,
                    'data': "owner",
                    'render': $.fn.dataTable.render.format_link()
                },
                {
                    'targets': 13,
                    'data': "permission",
                    'render': $.fn.dataTable.render.format_link()
                },
            ]
        });
        // When we've drawn the table, make the whole row clickable, links to href destination
        that._datatableObj.on("draw", function () {

            // create a fake element
            // add the rows HTML
            // Get the href from the Cell
            that._datatableObj.on("click", "tbody tr", function () {
                let el = $('<div></div>');
                el.html($(this)[0].outerHTML);

                window.location.href = $('a', el)[0].href
            });
        });
    }

}