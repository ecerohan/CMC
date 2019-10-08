function cmc_SDA= cmc(x,fs,P,levelnum)

%% This programm is used to extract constant-Q multi-level coefficients (CMC)
%% Jichen Yang, Rohan Kumar Das, Nina Zhou, "Extraction of octave spectra information for spoofing attack detection," 
%% IEEE/ACM Transactions on Audio, Speech and Language Processing, 2019.
%% Please cite to the above article while using this code.
%% Human Language Technology (HLT) Lab, Department of Electrical and Computer Engineering, National University of Singapore.
%% We would like to thank Massimiliano Todisco, Hector Delgado and Nicolas Evans who proposed constant-Q cepstral coefficient (CQCC) for anti-spoofing.
%% Our work is inspired by their work.
%% Please note that CQT_toolbox is required to run this code, which is available in the following link:
%% https://github.com/azraelkuan/asvspoof2017/tree/master/baseline/CQCC_v1.0/CQT_toolbox_2013
%% If you have some questions, please e-mail to us: NisonYoung@163.com, eleyji@nus.edu.sg, rohankd@nus.edu.sg. 
%% 2019-October-08.

%%Input:   x:        sampling points, 
%%         fs:       sampling rate, for ASVspoof 2015-2019, which equals 16KHz, 
%%         P:        the first top DCT result of every layer, in our experiments, we found that we can obtain the best performance when P equals 12.
%%         levelnum: number of levels (maximum value is 11)
%%Output:  CMC_SDA:  CMC feature with static, delta and acceleration, which includes static, delta and delta-delta. 

 
B = 96;                           %%  Frequency bin number in every octave
fmax = fs/2;                      %%  The maxmum frequency
oct = ceil(log2(fmax/20));        %%  Octave number set for CQT
fmin = fmax/2^oct;                %%  The minmimum frquency.
gamma = 228.7*(2^(1/B)-2^(-1/B)); %%  It is used in CQT, which can make the bandwidths equal a constant fraction of the ERB critical bandwidth


%%% CQT COMPUTATION
Xcq = cqt(x, B, fs, fmin, fmax, 'rasterize', 'full', 'gamma', gamma);    %% CQT

%%% LOG POWER SPECTRUM
absCQT = abs(Xcq.c);                                                     %% Magnitude spectrum
LogP_absCQT = log(absCQT.^2 + eps);                                      %% Octave power spectrum in log-scale

%%% CMC FEATURE EXTRACTION
oct1=dct(LogP_absCQT); oct1p=oct1(1:P,:);  oct1r=oct1(P+1:end,:);  %% the 1st layer tranform,  the 1st layer principal and residual infornation, p and r reprsent principal and residual
oct2=dct(oct1r);       oct2p=oct2(1:P,:);  oct2r=oct2(P+1:end,:);  %% The 2nd layer transform, the 2nd layer principal and residual information
oct3=dct(oct2r);       oct3p=oct3(1:P,:);  oct3r=oct3(P+1:end,:);
oct4=dct(oct3r);       oct4p=oct4(1:P,:);  oct4r=oct4(P+1:end,:);
oct5=dct(oct4r);       oct5p=oct5(1:P,:);  oct5r=oct5(P+1:end,:);
oct6=dct(oct5r);       oct6p=oct6(1:P,:);  oct6r=oct6(P+1:end,:);
oct7=dct(oct6r);       oct7p=oct7(1:P,:);  oct7r=oct7(P+1:end,:);
oct8=dct(oct7r);       oct8p=oct8(1:P,:);  oct8r=oct8(P+1:end,:);
oct9=dct(oct8r);       oct9p=oct9(1:P,:);  oct9r=oct9(P+1:end,:);
oct10=dct(oct9r);      oct10p=oct10(1:P,:);oct10r=oct10(P+1:end,:);
oct11=dct(oct10r);     oct11p=oct11(1:P,:);oct11r=oct11(P+1:end,:);%% The 11-th layer transform, the 11-th layer principal and residual information.

if(levelnum==1)
   CQcepstrum_temp=[oct1p];
elseif(levelnum==2)
   CQcepstrum_temp=[oct1p;oct2p];
elseif(levelnum==3)
   CQcepstrum_temp=[oct1p;oct2p;oct3p];
elseif(levelnum==4)
   CQcepstrum_temp=[oct1p;oct2p;oct3p;oct4p];
elseif(levelnum==5)
   CQcepstrum_temp=[oct1p;oct2p;oct3p;oct4p;oct5p];
elseif(levelnum==6)
   CQcepstrum_temp=[oct1p;oct2p;oct3p;oct4p;oct5p;oct6p];
elseif(levelnum==7)
   CQcepstrum_temp=[oct1p;oct2p;oct3p;oct4p;oct5p;oct6p;oct7p];
elseif(levelnum==8)
   CQcepstrum_temp=[oct1p;oct2p;oct3p;oct4p;oct5p;oct6p;oct7p;oct8p];
elseif(levelnum==9)
   CQcepstrum_temp=[oct1p;oct2p;oct3p;oct4p;oct5p;oct6p;oct7p;oct8p;oct9p];
elseif(levelnum==10)
   CQcepstrum_temp=[oct1p;oct2p;oct3p;oct4p;oct5p;oct6p;oct7p;oct8p;oct9p;oct10p];
else
   CQcepstrum_temp=[oct1p;oct2p;oct3p;oct4p;oct5p;oct6p;oct7p;oct8p;oct9p;oct10p;oct11p];
end
f_d = 1; % delta window size
fea = [CQcepstrum_temp; Deltas(CQcepstrum_temp,f_d); Deltas(Deltas(CQcepstrum_temp,f_d),f_d)]; %%Static+Delta+Acceleration
cmc_SDA=fea';
end


function D = Deltas(x,hlen)
% Delta and acceleration coefficients
% Reference:
%   Young S.J., Evermann G., Gales M.J.F., Kershaw D., Liu X., Moore G., Odell J., Ollason D.,
%   Povey D., Valtchev V. and Woodland P., The HTK Book (for HTK Version 3.4) December 2006.

win = hlen:-1:-hlen;
xx = [repmat(x(:,1),1,hlen),x,repmat(x(:,end),1,hlen)];
D = filter(win, 1, xx, [], 2);
D = D(:,hlen+1:(end - hlen));
D = D./(2*sum((1:hlen).^2));
end

