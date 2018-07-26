from django import template

register = template.Library()


@register.inclusion_tag('visualisation.html')
def render_vis():
    life = "life"
    return {"life": life}
