from django.shortcuts import render

def index(request):
    return render(request, 'web/index.html')

def current(request):
    return render(request, 'web/current_run.html')

def log_in(request):
    return render(request, 'web/log_in.html')

def private_index(request):
    return render(request, 'web/private_index.html')