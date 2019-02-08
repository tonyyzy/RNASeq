from django.shortcuts import render, redirect
from django.http import HttpResponse
from .forms import userRegistrationForm, userLoginForm
from django.contrib import messages # message.debug/info/success/warning/error
from django.contrib.auth import login, logout # imports the login funciton


def regiser_view(request):
    # return HttpResponse(" returned html directly.")
    if request.method == 'POST':
        form = userRegistrationForm(request.POST)
        if form.is_valid():
            print('\n FORM WAS VALID \n')
            user = form.save()
            login(request, user)
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


def login_view(request):
    # return HttpResponse('simple html')
    if request.method == 'POST':
        form = userLoginForm(data=request.POST)
        if form.is_valid():
            user = form.get_user()
            login(request, user)
            return redirect('analysis:samples_list')
    else:
        form = userLoginForm()
    return render(request, 'users/login.html', {'form': form})


def logout_view(request):
    # return HttpResponse('simple html')
    if request.method == 'POST':
        logout(request)
        return redirect('analysis:samples_list')
    else:
        form = userLoginForm()
    return render(request, 'users/login.html', {'form': form})
