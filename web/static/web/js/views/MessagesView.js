class MessagesView {

    constructor(htmlElement) {

        this._htmlElement = htmlElement;
    }

    update(model) {

            this._htmlElement.innerHTML = this._template(model);
    }

    _template(model) {
        
        return `
            <a href="#" class="dropdown-toggle" data-toggle="dropdown">
                <i class="fa fa-table"></i>
                <span class="label label-success">
                    <div id="all-messages-number">${model.length}</div>
                </span>
            </a>

            <ul class="dropdown-menu">
    
                <li>
                    <!-- inner menu: contains the messages -->
                    <ul class="menu">

                    ${model.map(n => {
                    return `
                            <li><a><span class="message">${DateHelper.datetimeToText(n.created_date)} - ${n.title}</span></a></li>
                        `
                    }).join('')}

                    </ul><!-- /.menu -->
                </li>
            </ul>
        `;
    }
}