from django.shortcuts import render, redirect
from django.http import HttpResponse
from .forms import userRegistrationForm
from django.contrib import messages # message.debug/info/success/warning/error

def regiser_view(request):
    # return HttpResponse(" returned html from users/views.py directly.")
    if request.method == 'POST':
        form = userRegistrationForm(request.POST)
        if form.is_valid():
            form.save()
            username = form.cleaned_data.get('username')
            messages.success(request, f'successfuly created: {username}')
            # args = {'username_link': username}
            # return render(request, 'register.html', args)
    else:
        form = userRegistrationForm()
    args = {'form': form}
    return render(request, 'register.html', args)
