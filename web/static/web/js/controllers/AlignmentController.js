/*global AlignmentController */
class AlignmentController {
    /**
     *
     * @param flowcellId {number} Primary key of flowcell record in database.
     */
    constructor(flowcellId) {
        this._flowcellId = flowcellId;
        this.coverageChartController = new CoverageChartController('coverage_div', "");
        this._axiosInstance = axios.create({
            headers: {"X-CSRFToken": getCookie('csrftoken')}
        });
        this.requestMappedChromosomes(flowcellId);
    }

    requestMappedChromosomes(flowcell_id) {
        /**
         * Request the chromosomes that have reads mapped to using minimap2
         * and update the select box on tab Mapping
         * TODO THIS NEEDS STREAMLINING
         */
        this._axiosInstance.get(`/api/v1/flowcells/${flowcell_id}/mapped-references`).then(response => {
            console.log(response.data)
            var select = document.getElementById('chromosome-id-select');
            var selected_index = select.selectedIndex;
            var selected_option = select[selected_index];

            while (select.length > 0) {
                select.remove(0);
            }

            var option = document.createElement('option');
            option.text = '--- Select ---';
            option.value = '-1';
            select.add(option);

            for (var i = 0; i < data.length; i++) {

                var option = document.createElement('option');

                var value_combination = data[i]['task_id'] + '_' + data[i]['barcode_name'] + '_' + data[i]['read_type_id'] + '_' + data[i]['chromosome_id'];

                option.text = data[i]['barcode_name'] + ' - ' + data[i]['read_type_name'] + ' - ' + data[i]['reference_name'] + ' - ' + data[i]['chromosome_name'];
                option.value = value_combination;

                if (selected_option === undefined) {

                } else {

                    if (option.text === selected_option.text) {
                        option.selected = 'selected';
                    }

                }

                select.add(option);
            }

        });

    }
}