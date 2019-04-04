from django import template

register = template.Library()

@register.filter(name='readable')
def readable(arg):
    arg = str(arg)
    abs_file = arg.split('/')
    readable_file = abs_file[-1]
    return readable_file

@register.filter(name='doubleMe')
def doubleMe(arg):
    return arg*2
