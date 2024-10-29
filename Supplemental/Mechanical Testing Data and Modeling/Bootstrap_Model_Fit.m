function [fits,stats] = Bootstrap_Model_Fit(lam,sig,indicator,to_plot)
%{
    Here model fits are bootstrapped to find the variability in the model
    parameters. lam and sig contain all of the stretch and engineering
    stress data (respectively), and indicator marks each of the data points
    with a integer value unique to the sample (positive values indicate
    compression tests and negative indicate tensile tests). To perform the
    bootstrap, random pairs of tensile and compression tests are paired
    together and the Yeoh model is fit to the combine compression and
    tensile data. This constitutes a single bootstrap sample of the model
    parameters.
%}

% get the IDs for the compression and tensile tests
compression_tests = unique(indicator(indicator > 0));
tensile_tests = unique(indicator(indicator < 0));

% Model Definition

% uniaxial deformation ( lam1 = lam, lam2 = lam3 = 1/sqrt(lam) )
I1 = @(lam) lam.^2 + 2./lam;

% Yeoh model Engineering stress
sig_Yeoh = @(lam,C1,C2,C3) 2*(lam - lam.^-2).*( C1 + 2*C2*(I1(lam) - 3) + 3*C3*(I1(lam)-3).^2);

N = 50; % how many compression/tensile combinations to create

C1s = zeros(N,1);
C2s = zeros(N,1);
C3s = zeros(N,1);

comp_test = datasample(compression_tests,N);
tens_test = datasample(tensile_tests,N);

% preallocating for saving all model fits
lam_all = linspace(min(lam),max(lam),50);
sig_all = zeros(N,50);

if to_plot
    figure;
end
for i=1:N
    LAM = [lam(indicator==comp_test(i)), lam(indicator==tens_test(i))];
    SIG = [sig(indicator==comp_test(i)), sig(indicator==tens_test(i))];

    if to_plot
        subplot(round(sqrt(N)),round(sqrt(N))+1,i); hold all
        plot(LAM,SIG,'.','Color',rand(3,1))
    end

    err = @(X) sum( (SIG - sig_Yeoh(LAM,X(1),X(2),X(3))).^2 ) + (1e10)*((X(1) < 0) + (X(3) < 0));
    X_ = fminsearch(err,[1,0,0]);
    C1 = X_(1);
    C2 = X_(2);
    C3 = X_(3);

    sig_all(i,:) = sig_Yeoh(lam_all,C1,C2,C3); 
    if to_plot
        plot(lam_all,sig_all(i,:),'--k')
    end

    C1s(i) = C1;
    C2s(i) = C2;
    C3s(i) = C3;

end

mean_sig = mean(sig_all,1);
std_sig = std(sig_all,[],1);

[C1_mean,C1_std] = RoundToSigFigs(mean(C1s),std(C1s));
[C2_mean,C2_std] = RoundToSigFigs(mean(C2s),std(C2s));
[C3_mean,C3_std] = RoundToSigFigs(mean(C3s),std(C3s));

% fprintf("Model Parameters: \n\tC1 = %.8f +/- %.8f,\n\tC2 = %.8f +/- %.8f, \n\tC3 = %.8f +/- %.8f\n",mean(C1s),std(C1s),mean(C2s),std(C2s),mean(C3s),std(C3s))
% fprintf("Rounded Model Parameters: \n\tC1 = %s +/- %s,\n\tC2 = %s +/- %s, \n\tC3 = %s +/- %s\n",num2str(C1_mean),num2str(C1_std),num2str(C2_mean),num2str(C2_std),num2str(C3_mean),num2str(C3_std))



if to_plot
    figure; hold all
    fill([lam_all,fliplr(lam_all)],[mean_sig + std_sig, fliplr(mean_sig - std_sig)],[0.8,0.2,0.6]*0.6,'FaceAlpha',0.2,'EdgeColor','none');
    plot(lam_all,mean_sig,'--','color',[0.8,0.2,0.6])
    scatter(lam,sig,2,'k','filled',"o","MarkerEdgeColor","none","MarkerFaceAlpha",0.7)
    xlabel("Stretch [mm/mm]")
    ylabel("Engineering Stress [MPa]")
    title("Yeoh Model Fit")
end

% calculate the R^2 for the mean fit to all of the data
sig_all_mod = sig_Yeoh(lam,mean(C1s),mean(C2s),mean(C3s));
SS_res = sum( ( sig_all_mod - sig ).^2 ); % sum of square residuals
SS_tot = sum( ( sig - mean(sig) ).^2 ); % total sum of squares
R2 = 1 - SS_res/SS_tot;


fits = struct();
fits.lam = lam;
fits.sig = sig;
fits.indicator = indicator;
fits.comp_samples = comp_test;
fits.tens_samples = tens_test;
fits.C1 = C1s;
fits.C2 = C2s;
fits.C3 = C3s;
fits.N = N;

stats = struct();
stats.C1_mean = mean(C1s);
stats.C1_std = std(C1s);
stats.C2_mean = mean(C2s);
stats.C2_std = std(C2s);
stats.C3_mean = mean(C3s);
stats.C3_std = std(C3s);
stats.lam_all = lam_all;
stats.sig_mean = mean_sig;
stats.sig_std = std_sig;
stats.R2 = R2;

end