function PD_data = loadNIRS_glucose_data(subject,eventT,eventN,storage)
    %PD_data = 1; 
    %subject
    %eventT
    %eventN
    %load('PD10_m_hypo_1.nirs','-mat');

    %load Glucose Data
    glucose_data = load(storage+":\PadovaPostDoc\BabyGluCo\formatted_PD_GuyPerkins\"+subject+".mat");


    %Load PD NIRS Data
    PD_data= load(storage+":\PadovaPostDoc\BabyGluCo\NIRS_data\"+subject+"\formatted\"+subject+"_"+eventT+"_"+eventN+".nirs",'-mat');
    
    %Combine PD Nirs Data and Glucose Data
    %Need to write NIRS data into PD data to make it general.
    PD_data.glucose = glucose_data;
    PD_data.SD =  getfield(PD_data,'PD_data',subject+"_"+eventT+"_X",'SD');
    PD_data.d =  getfield(PD_data,'PD_data',subject+"_"+eventT+"_X",'d');
    PD_data.t =  getfield(PD_data,'PD_data',subject+"_"+eventT+"_X",'t');
    PD_data.s =  getfield(PD_data,'PD_data',subject+"_"+eventT+"_X",'s');
    PD_data.aux =  getfield(PD_data,'PD_data',subject+"_"+eventT+"_X",'aux');
    %Need to write Glucose data into PD data to make it general.
    PD_data.glucose = getfield(PD_data,'glucose',subject);
end

%getfield(PD_data,'PD_data','PD10_m_hypo_X','SD')

%subject
%PD_data.PD_data.+subject+_m_hypo_X.SD

%subject+"_"+eventT+"_X"
%eventT
%eventN

%getfield(PD_data,'PD_data',subject+"_"+eventT+"_X",'SD')

%getfield(PD_data,'glucose','PD10')  