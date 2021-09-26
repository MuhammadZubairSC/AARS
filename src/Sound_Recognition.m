%Author:Muhammad Zubair SC
%Title: Automatic Sound Reconigtion System

function varargout = Sound_Recognition(varargin)
%SOUND_RECOGNITION MATLAB code for Sound_Recognition.fig
%
%      SOUND_RECOGNITION, by itself, creates a new SOUND_RECOGNITION or raises the existing
%      singleton*.
%
%      H = SOUND_RECOGNITION returns the handle to a new SOUND_RECOGNITION or the handle to
%      the existing singleton*.
%
%      SOUND_RECOGNITION('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in SOUND_RECOGNITION.M with the given input arguments.
%
%      SOUND_RECOGNITION('Property','Value',...) creates a new SOUND_RECOGNITION or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before Sound_Recognition_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to Sound_Recognition_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help Sound_Recognition

% Last Modified by GUIDE v2.5 03-Nov-2016 17:41:49

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @Sound_Recognition_OpeningFcn, ...
                   'gui_OutputFcn',  @Sound_Recognition_OutputFcn, ...
                   'gui_LayoutFcn',  [] , ...
                   'gui_Callback',   []);
if nargin && ischar(varargin{1})
    gui_State.gui_Callback = str2func(varargin{1});
end

if nargout
    [varargout{1:nargout}] = gui_mainfcn(gui_State, varargin{:});
else
    gui_mainfcn(gui_State, varargin{:});
end
% End initialization code - DO NOT EDIT


% --- Executes just before Sound_Recognition is made visible.
function Sound_Recognition_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to Sound_Recognition (see VARARGIN)
% Choose default command line output for Sound_Recognition
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes Sound_Recognition wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = Sound_Recognition_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


% --- Executes on selection change in pop_time_pre.
function pop_time_pre_Callback(hObject, eventdata, handles)
% hObject    handle to pop_time_pre (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns pop_time_pre contents as cell array
%        contents{get(hObject,'Value')} returns selected item from pop_time_pre
%Extracts the selected value from the uicontrol popup menu and store it in
%variable 'plot_1' and store it in container 'UserData', for out of the function
%use
plot_1=get(handles.pop_time_pre,'Value');
set(handles.pop_time_pre,'UserData',plot_1);

% --- Executes during object creation, after setting all properties.
function pop_time_pre_CreateFcn(hObject, eventdata, handles)
% hObject    handle to pop_time_pre (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in pop_wind_fre.
function pop_wind_fre_Callback(hObject, eventdata, handles)
% hObject    handle to pop_wind_fre (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns pop_wind_fre contents as cell array
%        contents{get(hObject,'Value')} returns selected item from pop_wind_fre
%Extracts the selected value from the uicontrol popup menu and store it in
%variable 'plot_2' and store it in container 'UserData', for out of the function
%use
plot_2=get(handles.pop_wind_fre,'Value');
set(handles.pop_wind_fre,'UserData',plot_2);


% --- Executes during object creation, after setting all properties.
function pop_wind_fre_CreateFcn(hObject, eventdata, handles)
% hObject    handle to pop_wind_fre (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in pop_Mel.
function pop_Mel_Callback(hObject, eventdata, handles)
% hObject    handle to pop_Mel (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns pop_Mel contents as cell array
%        contents{get(hObject,'Value')} returns selected item from pop_Mel
%Extracts the selected value from the uicontrol popup menu and store it in
%variable 'plot_3' and store it in container 'UserData', for out of the function
%use
plot_3=get(handles.pop_Mel,'Value');
set(handles.pop_Mel,'UserData',plot_3);

% --- Executes during object creation, after setting all properties.
function pop_Mel_CreateFcn(hObject, eventdata, handles)
% hObject    handle to pop_Mel (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in pop_mean.
function pop_mean_Callback(hObject, eventdata, handles)
% hObject    handle to pop_mean (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns pop_mean contents as cell array
%        contents{get(hObject,'Value')} returns selected item from pop_mean
%Extracts the selected value from the uicontrol popup menu and store it in
%variable 'plot_4' and store it in container 'UserData', for out of the function
%use
plot_4=get(handles.pop_mean,'Value');
set(handles.pop_mean,'UserData',plot_4);

% --- Executes during object creation, after setting all properties.
function pop_mean_CreateFcn(hObject, eventdata, handles)
% hObject    handle to pop_mean (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in quit.
function quit_Callback(hObject, eventdata, handles)
% hObject    handle to quit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
%Wait for the user to give input
%closes all the open windows
uiwait(msgbox('Thank you for using Automatic Aircraft Recognition System','Exit','custom',imread('Exit.png'))); %Ending statement of program
pause(5);
close(gcbf);
close('all'); 

% --- Executes on button press in test.
function test_Callback(hObject, eventdata, handles)
% hObject    handle to test (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
%% Initialization section
%Refreshes all the plot in the GUI
cla(handles.axes1); 
cla(handles.axes2); 
cla(handles.axes3); 
cla(handles.axes4);
%Retrieves the value of different parameters from uicontrols of Automatic
%Speech Recognition System
plot_1=get(handles.pop_time_pre,'UserData');
plot_2=get(handles.pop_wind_fre,'UserData');
plot_3=get(handles.pop_Mel,'UserData');
plot_4=get(handles.pop_mean,'UserData');
a_y=get(handles.analysis,'Value');
d_t=get(handles.detection,'Value');
audio_signal_selection= get(handles.audio_signals,'SelectedObject');
set(handles.text_display,'String','');
%% Error Handling Section
%If the user chooses Analysis or both detection and analysis fo sound to be
%done
if(a_y==1 ||(a_y==1 && d_t==1))
%If the variable is not a number, an alphabet perhaps or contains no value
%then error message will pop up in order to guide the user back to the main
%screen to define appropriate selection from each popup menu
if((isempty(plot_1)) || (isnan(plot_1))||(isempty(plot_2)) || (isnan(plot_2))||(isempty(plot_3))|| (isnan(plot_3))||(isempty(plot_4)) || (isnan(plot_4)))
    msgbox('Please select any one from the given options for each graph popup menu.','Error Message','custom',imread('errorsign.png'));
    ploterror=audioread('plots_error.mp3');sound(ploterror,44100/2);
    set(handles.test,'Enable','off');
    set(handles.train,'Enable','off'); 
    set(handles.reset,'Enable','off'); 
    set(handles.quit,'Enable','off'); 
    pause(3);
    set(handles.test,'Enable','on');
    set(handles.train,'Enable','on'); 
    set(handles.reset,'Enable','on'); 
    set(handles.quit,'Enable','on'); 
    return;
end
end
%If user does not define the function of GUI at all the follwing program
%rectify it
if(d_t==0 && a_y==0)
 msgbox('Please specify the function of GUI','Error Message','custom',imread('errorsign.png'));
 return;
end
%Checking  whether the excel spreadsheet has data in it if not the program
%ask user to tain the system first
if((isempty(xlsread('Train.xls','Sheet4','A2:K342'))))
    msgbox('Please train the system before testing it.','Error Message','custom',imread('errorsign.png'));
    trainbef=audioread('Train before Test.mp3');
    sound(trainbef,44100/2);
    %Turn off all the buttons and popup menus
    set(handles.test,'Enable','off');
    set(handles.train,'Enable','off'); 
    set(handles.reset,'Enable','off'); 
    set(handles.quit,'Enable','off'); 
    pause(5);
    %Turn On all the buttons and popup menus
    set(handles.test,'Enable','on');
    set(handles.train,'Enable','on'); 
    set(handles.reset,'Enable','on'); 
    set(handles.quit,'Enable','on'); 
    return;
end
% Turning off all the buttons and popup menus so that user does not change
% any parameter during simulation
set(handles.test,'Enable','off');
set(handles.train,'Enable','off'); 
set(handles.reset,'Enable','off'); 
set(handles.quit,'Enable','off');
set(handles.pop_mean,'Enable','off');
set(handles.pop_Mel,'Enable','off'); 
set(handles.pop_time_pre,'Enable','off'); 
set(handles.pop_wind_fre,'Enable','off');
%Extracting the string from radio panel variable 'audio_signal_selection'
audio_S=get(audio_signal_selection,'String');
%Taking random values and put in two different variables to provide
%probability for getting different samples of same signal
shift_a=rand(1,1);
shift_b=rand(1,1);
%Switch statement for selecting the audio signal based on the given data
switch(audio_S)
     %In this case of 'shift_a' is greater than 'shift_b' then the original
     %signal will be used and if its not then the sampled signal will be
     %used as the input signal
    case 'B-25 Mitchell Jet'
     if(shift_a>shift_b)
     [Input,fs]=audioread('B-25 Mitchell.wav');
     else
     [Input,fs]=audioread('B-25 Mitchell (sample 1).wav');
     end
      %In this case of 'shift_a' is greater than 'shift_b' then the original
     %signal will be used and if its not then the sampled signal will be
     %used as the input signal
    case 'Cessna Private Jet'
     if(shift_a>shift_b)
     [Input,fs]=audioread('Cessna.wav');
     else
     [Input,fs]=audioread('Cessna (sample 1).wav');
     end
    
     %In this case of 'shift_a' is greater than 'shift_b' then the original
     %signal will be used and if its not then the sampled signal will be
     %used as the input signal
     
    case 'F-16 Fighter Falcon'
     if(shift_a>shift_b)
     [Input,fs]=audioread('F-16 Fighting Falcon.wav');
     else
     [Input,fs]=audioread('F-16 fighting falcon (sample 1).wav');
     end
     
     %In this case of 'shift_a' is greater than 'shift_b' then the original
     %signal will be used and if its not then the sampled signal will be
     %used as the input signal
     
     case 'Helicopter'
     if(shift_a>shift_b)
     [Input,fs]=audioread('Helicopter.wav');
     else
     [Input,fs]=audioread('Helicopter (sample 1).mp3');
     end
     
     %In this case of 'shift_a' is greater than 'shift_b' then the original
     %signal will be used and if its not then the sampled signal will be
     %used as the input signal
     %Playing the sound of selected signal
     % wait for 7 seconds
    case 'MIG-21 Fighting Jet'
     if(shift_a>shift_b)
     [Input,fs]=audioread('Mig-21.wav');
     else
     [Input,fs]=audioread('Mig-21 (sample 1).wav');
     end
     
    otherwise
     %The case takes one arbitrary signal of Jet and uses it as an input if
     %none of the above cases exists
     %Playing the sound of selected signal
     % wait for 10 seconds
     [Input,fs]=audioread('FA-18 Fighter Jet.wav');
        
end
%Extract the mono section of sterero sound assign it to variable 'audr'
%Takes the length of it and takes the time vector
Ip=Input(:,1);
audr=Ip;
N=length(audr);
time=linspace(0,N/fs,N);
%The Transfer function of per-emphasis filter
%Creates the filter based on inserted Transfer function
%Takes frame duration and round the resulting frame length ot get the whole
%value and record the overlap period
%Initialize the matircies before modifying them in loop
%Take the neares t power of two which compensates framelength and creates
%the Fourier transfrom vairbale having that length equal to 2^nfft.
%assign 26 filterbanks to variable F
%Initialization of Mel Filter bank variable and creating hamming window
%function
preEmph=2*[1 -1];
aud = filter(preEmph,1,audr);
frameduration=0.025;
framelen=round(frameduration*fs);
overlap=round(0.01*fs);
lasframe=floor(N/framelen);
audf=zeros(framelen,lasframe);
filterd=zeros(framelen,lasframe);
nfft=nextpow2(framelen);
XFTful=zeros(2^nfft,lasframe);
fft_length=2^nfft/2;
Xft=zeros((2^nfft)/2,lasframe);
F=26;
M=zeros(F,(2^nfft)/2);
k=0:(2^nfft)/2-1;
wn=hamming(framelen);
%The loop divides the audio in terms of frames of calculated frame length
%having same overlapping as calculated and apply the hamming window ot each
%frame and then takes the FFT of each frame and extracts the magnitude of it and discrad the half of FFT
%spectrum
 for o=1:lasframe
 audf(:,o)=aud(((o-1)*(framelen-overlap)+1):(o*framelen-(o-1)*overlap));
 filterd(:,o)=wn.*audf(:,o);
 XFTful(:,o)=abs(fft(filterd(:,o),2^nfft)); %for nfft
 Xft(:,o)=XFTful(1:(2^nfft)/2,o);
 end
 %Calculates the overall FFT of the audio signal and takes the magnitude of
 %it and extracts the half of the spectrum and creates the vector shwoning
 %the frequency in Hertz 
 nfft2=2^(nextpow2(length(aud)));
 OvXft=abs(fft(aud,nfft2));
 OvXft2=OvXft(1:nfft2/2);
 fd=(0:N-1).*fs/nfft2;
 
 lfk=k*fs/framelen;%linear frequncy
 lfmax=fs;%Maximum feqeuncy should be less than nuqyist rate fs/2; 
 lfmin=100;%miminum suggestable frequency in >100 
 phimax=2595*log10(1+lfmax/700);%maximum mel frequnecy
 phimin=2595*log10(1+lfmin/700);%minimum mel frequnecy
 Dphif=(phimax-phimin)/(F+1);%fixed freqeuncy resolution in mel frequency scale
 phifC=(1:F)*Dphif;%Creating a vector for center frequency in MEL 
 lfc=700*(10.^(phifC/2595)-1);%Converts the vector to time 
 %Implementation of Mel Filter Bank formula
 %The first loop result in 26 rows and seconds one result in column  equal to half of length of FFT spectrum matrix 
 for m=2:F-1;
     for k=1:(2^nfft)/2
 if((lfk(k)<lfc(m-1))||(lfk(k)>=lfc(m+1)))
     M(m,k)=0;
elseif ((lfc(m-1)<=lfk(k))&&(lfk(k)<lfc(m)))
M(m,k)=(lfk(k)-lfc(m-1))/(lfc(m)-lfc(m-1));
elseif ((lfc(m)<=lfk(k))&&(lfk(k)<lfc(m+1)))
M(m,k)=(lfk(k)-lfc(m+1))/(lfc(m)-lfc(m+1));
end
    end
 end
%Initialzing the matrix before using them in loop
EHalf=zeros((2^nfft)/2,1);
Energy=zeros(F,lasframe);
%Mel frequncy wraping (multiplying the FFT of each frame with Mel filter
%First loop calculates the coeccfients for each frame and second loop
%increases the frame number to the last frame
for o=1:lasframe
    for m=1:F;
    for n=1:(2^nfft)/2
 EHalf(n)=M(m,n)*Xft(n,o);
    end
    Energy(m,o)=sum(EHalf);
    EHalf=zeros((2^nfft)/2,1);
    end
end
%Takes the logarithm of the naswer to get lof MFCC
logEnergy=log(Energy);
%DCT program
%Countering the inf value due to filteration by comparing each element of
%Log MFCC matrix
for o=1:lasframe
    for m=1:26
     if((logEnergy(m,o)==Inf)||(logEnergy(m,o)==-Inf)||isnan(logEnergy(m,o)))
         logEnergy(m,o)=0;
     end
    end
end
%Takes the DCT (Discrete cosine transfrom) of Log MFCC
%Sum the entire columns to get only one column vector of mean MFCC and
%divdes the result by total number of frames
cn=dct(logEnergy);
cnoverall=sum(cn');
cnoverall=cnoverall';
MFCC_Ip=round(cnoverall/lasframe);
%Takes the log Nenergy of signal in time and counter any inf value by
%noting the index of it and flusing them off from the matrix
logE_Input=log(sum(audf.^2))';
Inf_T=isinf(logE_Input);
logE_Input(find(Inf_T))=0;
%If analysis is required then only the plots will be shown
if(a_y==1)
%Switch structure for displaying the user selected graph on ‘axes1'
%the 1st case is ceased, 2nd, 3rd, 4th and 5th shows the
%Time Domain Representation of Input Signal , Pole-Zero Plot of Pre-emphasis Filter,Magnitude Response of Pre-emphasis Filter
%and Transfer function of Pre-emphasis Filter respectively
switch(plot_1)
    case 2
        axes(handles.axes1);
        plot(time,audr);
        xlabel('Time in (Seconds)'),ylabel('Amplitude'),title('Time Domain Representation of Input Signal');
    case 3
        axes(handles.axes1);
        zplane(preEmph,(1)); 
        title('Pole-Zero Plot of Pre-emphasis Filter');
    case 4
        axes(handles.axes1);
        [h,w]=freqz(preEmph,(1));
        plot(w/pi,20*log10(abs(h)));
        title('Magnitude Response of Pre-emphasis Filter'),xlabel('Normalized Frequency (\times\pi rad/sample)'),ylabel('Magnitude (dB)');
    case 5
        axes(handles.axes1);
        tfestimate(audr,aud),title('Transfer function of Pre-emphasis Filter');
    otherwise
end

%Switch structure for displaying the user selected graph on ‘axes2'
%the 1st case is ceased, 2nd, 3rd, 4th and 5th shows the
%Magnitude Response of Hamming Window , Overall frequency spectrum of the signal,Framed Audio signal
%and Framed Frequnecy Spectrum of Input signal respectively by calculating
%their x axis vector similar to pervious calculation
switch(plot_2)
    case 2
        axes(handles.axes2);
        [h1,w1]=freqz(wn);
        plot(w1/pi,20*log10(abs(h1)));
        title('Magnitude Response of Hamming Window'),xlabel('Normalized Frequency (\times\pi rad/sample)'),ylabel('Magnitude (dB)');
    case 3
        axes(handles.axes2);
        plot(fd(1:nfft2/2),OvXft2);
        title('Overall frequency spectrum of the signal'),xlabel('Frequency in (Hertz)'),ylabel('Magnitude');
    case 4
        axes(handles.axes2);
        plot(audf);
        title('Framed Audio signal'),xlabel('Frame length'),ylabel('Amplitude');
    case 5
        axes(handles.axes2);
        frequency_axis=zeros((2^nfft)/2,lasframe);
        for l=1:lasframe
            frequency_axis(:,l)=linspace(0,fs/2,(2^nfft)/2);
        end
        plot(frequency_axis,Xft);
        title('Framed Frequnecy Spectrum of Input signal'),xlabel('Frequnecy in (Hertz)'),ylabel('Magnitude');
    otherwise
end
%Switch structure for displaying the user selected graph on ‘axes3'
%the 1st case is ceased, 2nd, 3rd, 4th and 5th shows the
%Mel Filter Bank Plot , Logarithmic Mel Frequncy Cepstral Coefficient Plot,Mel Frequncy Cepstral Coefficient (MFCC) Time domain Plot
%and Fransfer Function of Mel Filter Bank respectively
switch(plot_3)
    case 2
        axes(handles.axes3);
        plot(linspace(0,fs/2,length(M)),M);
        title('Mel Filter Bank Plot'),xlabel('Frequency in (Hertz)'),ylabel('Amplitude');
    case 3
        axes(handles.axes3);
        plot(linspace(0,fs/2,F),logEnergy);
        title('Logarithmic Mel Frequncy Cepstral Coefficient Plot'),xlabel('Frequency in (Hertz))'),ylabel('Magnitude(dB)');
    case 4
        axes(handles.axes3);
        plot(linspace(0,N/fs,F),cn);
        title('Mel Frequncy Cepstral Coefficient (MFCC) Time domain Plot'),xlabel('Time in (Seconds))'),ylabel('MFCC');
    case 5
        axes(handles.axes3);
        Energy_pad=Energy;
        Energy_pad(F+1:fft_length,:)=zeros(fft_length-F,lasframe);%zero padding
        tfestimate(Xft,Energy_pad,fft_length,[],[],fs);
        title('Transfer Function of Mel Filter Bank');
    otherwise
end

%Switch structure for displaying the user selected graph on ‘axes4'
%the 1st case is ceased, 2nd, 3rd, 4th and 5th shows the
%Mean  Mel Frequncy Cepstral Coefficient Plot and Plot for Logarithmic Energy in Time domain ,respectively
switch(plot_4)
    case 2
        axes(handles.axes4);
        plot(linspace(0,N/fs,F),MFCC_Ip);
        title('Mean  Mel Frequncy Cepstral Coefficient Plot'),xlabel('Time in (Seconds)'),ylabel('Mean MFCC Amplitude');
   
    case 3
        axes(handles.axes4);
        plot(linspace(0,N/fs,length(logE_Input)),logE_Input);
        title('Plot for Logarithmic Energy in Time domain'),xlabel('Time in (Seconds)'),ylabel('Logarithmic Energy (dB)');
    otherwise
end
end

%Detection has been chosen and both analysis and detection has been
%chosen by the user
%Takes the round of log Energy
%Takes the data from the Excel file 'Train.xls' and extract each column
%from the data and assign it to the relevant airplane genre
if((a_y==1 && d_t==1)|| d_t==1)
logE_Input=round(logE_Input);
data= xlsread('Train.xls','Sheet4','A2:T419');
A_mitchell=round(data(:,1));
B_cessna=round(data(:,2));
C_falcon=round(data(:,3));
D_helicopter=round(data(:,4));
E_mig=round(data(:,5));
F_E_mitchell=round(data(:,6));
%Checks if there is any value which is not a number and gives the index of
%these elements, which thereby clips the element off
Nan_F=isnan(F_E_mitchell);
F_E_mitchell(find(Nan_F))=[];
G_E_cessna=round(data(:,7));
Nan_G=isnan(G_E_cessna);
G_E_cessna(find(Nan_G))=[];
H_E_falcon=round(data(:,8));
Nan_H=isnan(H_E_falcon);
H_E_falcon(find(Nan_H))=[];
I_E_helicopter=round(data(:,9));
Nan_I=isnan(I_E_helicopter);
I_E_helicopter(find(Nan_I))=[];
J_E_mig21=round(data(:,10));
Nan_J=isnan(J_E_mig21);
J_E_mig21(find(Nan_J))=[];
K_mitchell=round(data(:,11));
L_cessna=round(data(:,12));
M_falcon=round(data(:,13));
N_helicopter=round(data(:,14));
O_mig=round(data(:,15));
P_E_mitchell=round(data(:,16));
Nan_P=isnan(P_E_mitchell);
P_E_mitchell(find(Nan_P))=[];
Q_E_cessna=round(data(:,17));
Nan_Q=isnan(Q_E_cessna);
Q_E_cessna(find(Nan_Q))=[];
R_E_falcon=round(data(:,18));
Nan_R=isnan(R_E_falcon);
R_E_falcon(find(Nan_R))=[];
S_E_helicopter=round(data(:,19));
Nan_S=isnan(S_E_helicopter);
S_E_helicopter(find(Nan_S))=[];
T_E_mig21=round(data(:,20));
Nan_T=isnan(T_E_mig21);
T_E_mig21(find(Nan_T))=[];
%If the random variable comparison is positive then only the cross
%correlation will be done
if(shift_a>shift_b)
%Reads the audiofile and takes the mono part from it.
y1=audioread('B-25 Mitchell.wav');
y1=y1(:,1);
y2=audioread('Cessna.wav');
y2=y2(:,1);
y3=audioread('F-16 Fighting Falcon.wav');
y3=y3(:,1);
y4=audioread('Helicopter.wav');
y4=y4(:,1);
y5=audioread('Mig-21.wav');
y5=y5(:,1);
%Calculates the maximum value of corss correlation and assign it to the
%matrix 'Sample_xcor'
a1=max(xcorr(Ip,y1));
a2=max(xcorr(Ip,y2));
a3=max(xcorr(Ip,y3));
a4=max(xcorr(Ip,y4));
a5=max(xcorr(Ip,y5));
Sample_xcor=[a1;a2;a3;a4;a5];
%For loop for comparing each element of array with threshold value and then displaying the result for containing the code word for each jet detection 
for h=1:5;
    if(Sample_xcor(h)>=37000)
a=1;
    elseif((Sample_xcor(h)>=3230)&&(Sample_xcor(h)<=4000))
a=2;
    elseif((Sample_xcor(h)>=24800)&&(Sample_xcor(h)<=24900))
a=3;
    elseif(Sample_xcor(h)<=0.0004)
a=4;
    elseif((Sample_xcor(h)>=9840)&&(Sample_xcor(h)<=9850))
a=5;    
    end
end
%If the random variable comparison gives negative value then no corss
%correlation will be done
else
a=0;
end
%If random variable gives positive value
if(shift_a>shift_b)
%Compares the sum of vector obtained by comparing the input MFCC with the
%database MFCC and checks if it is greater than and equal to 12 together
%with the similar approach for log energy of input with database and the
%autocorrelation varaible comparison
%If any of the case equuals it shows that particular jet has been detected
%via narration effect
if(((sum(MFCC_Ip(1:13)==A_mitchell(1:13)))>=12)&&(a==1)&&(sum(logE_Input==F_E_mitchell)>=length(logE_Input)))
    set(handles.text_display,'String','B-25 Mitchell Fighter Jet has been detected');
    detec=audioread('B-25 Mitchell_detected.mp3');sound(detec,44100/2);
    pause(5);
elseif(((sum(MFCC_Ip(1:13)==B_cessna(1:13)))>=12)&&(a==2)&&(sum(logE_Input==G_E_cessna)>=length(logE_Input)))
    set(handles.text_display,'String','Cessna Private Jet has been detected');
    detec=audioread('Cessna_detected.mp3');sound(detec,44100/2);
    pause(3);
elseif(((sum(MFCC_Ip(1:13)==C_falcon(1:13)))>=12)&&(a==3)&&(sum(logE_Input==H_E_falcon)>=length(logE_Input)))
    set(handles.text_display,'String','F-16 Fighting Falcon Jet has been detected');
    detec=audioread('F-16_detected.mp3');sound(detec,44100/2);
    pause(3);
elseif(((sum(MFCC_Ip(1:13)==D_helicopter(1:13)))>=12)&&(a==4)&&(sum(logE_Input==I_E_helicopter)>=length(logE_Input)))
    set(handles.text_display,'String','Private Helicopter has been detected');
    detec=audioread('Helicopter_detected.mp3');sound(detec,44100/2);
    pause(3);
elseif(((sum(MFCC_Ip(1:13)>=E_mig(1:13)))>=12)&&(a==5)&&(sum(logE_Input==J_E_mig21)>=length(logE_Input)))
    set(handles.text_display,'String','MIG-21 Fighter Jet has been detected');
    detec=audioread('Mig-21_detected.mp3');sound(detec,44100/2);
    pause(3);
else
    set(handles.text_display,'String','None from the database');
end
else
%Compares the sum of vector obtained by comparing the input MFCC with the
%database MFCC and checks if it is greater than and equal to 12 together
%with the similar approach for log energy of input with database. Here
%cross correlation is not done.
  if(((sum(MFCC_Ip(1:13)==K_mitchell(1:13)))>=12)&&(a==0)&&(sum(logE_Input==P_E_mitchell)>=length(logE_Input)))
    set(handles.text_display,'String','B-25 Mitchell Fighter Jet has been detected');
    detec=audioread('B-25 Mitchell_detected.mp3');sound(detec,44100/2);
    pause(5);
elseif(((sum(MFCC_Ip(1:13)==L_cessna(1:13)))>=12)&&(a==0)&&(sum(logE_Input==Q_E_cessna)>=length(logE_Input)))
    set(handles.text_display,'String','Cessna Private Jet has been detected');
    detec=audioread('Cessna_detected.mp3');sound(detec,44100/2);
    pause(3);
elseif(((sum(MFCC_Ip(1:13)==M_falcon(1:13)))>=12)&&(a==0)&&(sum(logE_Input==R_E_falcon)>=length(logE_Input)))
    set(handles.text_display,'String','F-16 Fighting Falcon Jet has been detected');
    detec=audioread('F-16_detected.mp3');sound(detec,44100/2);
    pause(3);
elseif(((sum(MFCC_Ip(1:13)==N_helicopter(1:13)))>=12)&&(a==0)&&(sum(logE_Input==S_E_helicopter)>=length(logE_Input)))
    set(handles.text_display,'String','Private Helicopter has been detected');
    detec=audioread('Helicopter_detected.mp3');sound(detec,44100/2);
    pause(3);
elseif(((sum(MFCC_Ip(1:13)>=O_mig(1:13)))>=12)&&(a==0)&&(sum(logE_Input==T_E_mig21)>=length(logE_Input)))
    set(handles.text_display,'String','MIG-21 Fighter Jet has been detected');
    detec=audioread('Mig-21_detected.mp3');sound(detec,44100/2);
    pause(3);
else
    set(handles.text_display,'String','None from the database');
  end
end
end
%Tuns on all the buttons and popup menu after completing the simulation
set(handles.test,'Enable','on');
set(handles.train,'Enable','on'); 
set(handles.reset,'Enable','on'); 
set(handles.quit,'Enable','on');
set(handles.pop_mean,'Enable','on');
set(handles.pop_Mel,'Enable','on'); 
set(handles.pop_time_pre,'Enable','on'); 
set(handles.pop_wind_fre,'Enable','on');



% --- Executes on button press in reset.
function reset_Callback(hObject, eventdata, handles)
% hObject    handle to reset (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
%closes the current base figure
%opens the main simulation function again 
close(gcbf); 
pause(3);
Sound_Recognition();

% --- Executes on button press in pushbutton4.
function pushbutton4_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton4 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on button press in train.
function train_Callback(hObject, eventdata, handles)
% hObject    handle to train (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
%Clears the text on the screen and switch off all the buttons and the popup
%menus
set(handles.text_display,'String',' ');
set(handles.test,'Enable','off');
set(handles.train,'Enable','off'); 
set(handles.reset,'Enable','off'); 
set(handles.quit,'Enable','off');
set(handles.pop_mean,'Enable','off');
set(handles.pop_Mel,'Enable','off'); 
set(handles.pop_time_pre,'Enable','off'); 
set(handles.pop_wind_fre,'Enable','off');
%For loop for extracting the information from five signals, one at a time.
for h=0:4
    %Takes the audio signal and extract the mono from it and extracts the
    %Mean MFCC and Log energy features by using function
    %'mfcc_program_mean' and 'Energy_program' function and then wirting
    %them in file the same steps proceed for the sampled signal.
    if(h==0)
    [audA, fs]=audioread('B-25 Mitchell.wav');
    audr=(audA(:,1));
    cnoverall=mfcc_program_mean(audr, fs);
    logE=Energy_program(audr, fs)';
    xlswrite('Train.xls',logE,'Sheet4','F2:F342');
    xlswrite('Train.xls',cnoverall,'Sheet4','A2:A15');
    [audA, fs]=audioread('B-25 Mitchell (sample 1).wav');
    audr=(audA(:,1));
    cnoverall=mfcc_program_mean(audr, fs);
    logE=Energy_program(audr, fs)';
    xlswrite('Train.xls',logE,'Sheet4','P2:P419');
    xlswrite('Train.xls',cnoverall,'Sheet4','K2:K15');
    elseif(h==1)
    [audA, fs]=audioread('Cessna.wav');
    audp=(audA(:,1));
    cnoverall=mfcc_program_mean(audp, fs);
    logE=Energy_program(audp, fs)';
    xlswrite('Train.xls',logE,'Sheet4','G2:G55');
    xlswrite('Train.xls',cnoverall,'Sheet4','B2:B15');
    [audA, fs]=audioread('Cessna (sample 1).wav');
    audr=(audA(:,1));
    cnoverall=mfcc_program_mean(audr, fs);
    logE=Energy_program(audr, fs)';
    xlswrite('Train.xls',logE,'Sheet4','Q2:Q419');
    xlswrite('Train.xls',cnoverall,'Sheet4','L2:L15');
    elseif(h==2)
    [audA, fs]=audioread('F-16 Fighting Falcon.wav');
    audp=(audA(:,1));
    cnoverall=mfcc_program_mean(audp, fs);
    logE=Energy_program(audp, fs)';
    xlswrite('Train.xls',logE,'Sheet4','H2:H350');
    xlswrite('Train.xls',cnoverall,'Sheet4','C2:C15');
    [audA, fs]=audioread('F-16 fighting falcon (sample 1).wav');
    audr=(audA(:,1));
    cnoverall=mfcc_program_mean(audr, fs);
    logE=Energy_program(audr, fs)';
    xlswrite('Train.xls',logE,'Sheet4','R2:R419');
    xlswrite('Train.xls',cnoverall,'Sheet4','M2:M15');
    elseif(h==3)
    [audA, fs]=audioread('Helicopter.wav');
    audp=(audA(:,1));
    cnoverall=mfcc_program_mean(audp, fs);
    logE=Energy_program(audp, fs)';
    xlswrite('Train.xls',logE,'Sheet4','I2:I419');
    xlswrite('Train.xls',cnoverall,'Sheet4','D2:D15');
    [audA, fs]=audioread('Helicopter (sample 1).mp3');
    audr=(audA(:,1));
    cnoverall=mfcc_program_mean(audr, fs);
    logE=Energy_program(audr, fs)';
    xlswrite('Train.xls',logE,'Sheet4','S2:S419');
    xlswrite('Train.xls',cnoverall,'Sheet4','N2:N15');
    else
    [audA, fs]=audioread('Mig-21.wav');
    audp=(audA(:,1));
    cnoverall=mfcc_program_mean(audp, fs);
    logE=Energy_program(audp, fs)';
    xlswrite('Train.xls',logE,'Sheet4','J2:J268');
    xlswrite('Train.xls',cnoverall,'Sheet4','E2:E15');
    [audA, fs]=audioread('Mig-21 (sample 1).wav');
    audr=(audA(:,1));
    cnoverall=mfcc_program_mean(audr, fs);
    logE=Energy_program(audr, fs)';
    xlswrite('Train.xls',logE,'Sheet4','T2:T419');
    xlswrite('Train.xls',cnoverall,'Sheet4','O2:O15');
    end
end
%displays that the system has been trained and switches on all the buttons
%and popupmenus
set(handles.text_display,'String','The System has been trained, successfully.');
set(handles.test,'Enable','on');
set(handles.train,'Enable','on'); 
set(handles.reset,'Enable','on'); 
set(handles.quit,'Enable','on');
set(handles.pop_mean,'Enable','on');
set(handles.pop_Mel,'Enable','on'); 
set(handles.pop_time_pre,'Enable','on'); 
set(handles.pop_wind_fre,'Enable','on');


% --- Executes on button press in detection.
function detection_Callback(hObject, eventdata, handles)
% hObject    handle to detection (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of detection


% --- Executes on button press in analysis.
function analysis_Callback(hObject, eventdata, handles)
% hObject    handle to analysis (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of analysis
