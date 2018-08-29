from django import template

register = template.Library()


@register.inclusion_tag('visualisation.html')
def render_vis():
    """
    a function that returns the above template, used on the metagenomics.html
    page contains all the metagenomics visualisations
    :return:
    """
    life = "life"
    return {"life": life}
