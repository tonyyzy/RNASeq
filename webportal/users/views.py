from django.shortcuts import render, redirect
from django.http import HttpResponse
from .forms import userRegistrationForm
from django.contrib import messages # message.debug/info/success/warning/error

def regiser_view(request):
    # return HttpResponse(" returned html directly.")
    if request.method == 'POST':
        form = userRegistrationForm(request.POST)
        if form.is_valid():
            print('\n FORM WAS VALID \n')
            form.save()
            username = form.cleaned_data.get('username')
            messages.success(request, f'successfuly created: {username}')
            # args = {'username_link': username}
            # return render(request, 'register.html', args)
            # return redirect('analysis:upload_session')
            return redirect('users:regiser_view')

    else:
        print(f'\n getting form with method: {request.method} \n')
        form = userRegistrationForm()
    args = {'form': form}
    return render(request, 'users/signup.html', args)
