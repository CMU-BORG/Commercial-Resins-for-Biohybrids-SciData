%% Biocompatibility Mechanical Testing Model Fitting
% Last Updated: MB 5/13/2024

clear all; close all; clc

%% Load Data
saveDir = "MATLAB Outputs\20250212\";
masterDataFile = "masterRawData_clean.xlsx";

out_file = saveDir + "Material Properties.xlsx";

%{
Useful Links:
importOptions/Preview: https://www.mathworks.com/help/matlab/ref/matlab.io.text.delimitedtextimportoptions.preview.html
combining categorical arrays: https://www.mathworks.com/help/matlab/matlab_prog/combine-categorical-arrays-using-multiplication.html

%}
geomOpts = detectImportOptions(masterDataFile,"Sheet","SpecificGeometries");
geomOpts.VariableTypes(5:6) = {'string','string'};
compOpts = detectImportOptions(masterDataFile,"Sheet","Compression_Rem");
compOpts.VariableTypes(4:5) = {'string','string'};
tensOpts = detectImportOptions(masterDataFile,"Sheet","Tensile_Rem");
tensOpts.VariableTypes(4:5) = {'string','string'};

geom_data = readtable(masterDataFile,geomOpts);
compress_data = readtable(masterDataFile,compOpts);
tensile_data = readtable(masterDataFile,tensOpts);

geom_data.Resin = lower(geom_data.Resin);
compress_data.Resin = lower(compress_data.Resin);
tensile_data.Resin = lower(tensile_data.Resin);

resins = unique(geom_data.Resin);
tests = unique(geom_data.TestType);
sterilization = {'NS','Autoclave','EtOH'}; %unique(geom_data.Sterilization);

N_interp_points = 50;
markers = ["o","s","^",">","<","v","*","square","diamond","+","pentagram","hexagram"];


%% Stress-Stretch Model

% stretch when accounting for tare load displacement
LAM = @(dL, a, L0) 1 + (dL + a)/L0;

% uniaxial deformation ( lam1 = lam, lam2 = lam3 = 1/sqrt(lam) )
I1 = @(lam) lam.^2 + 2./lam;

% Yeoh model Engineering stress
sig_Yeoh = @(lam,C1,C2,C3) 2*(lam - lam.^-2).*( C1 + 2*C2*(I1(lam) - 3) + 3*C3*(I1(lam)-3).^2);

% colors = {1.25*[0.2,0.8,0.6], [0.95,0.33,0.64], [0.09,0.45,0.6]};
colors = {[80,215,50]/255,[253,46,122]/255,[78,95,125]/255};
alpha = 0.3;
lines = {'-','--','-.'};
markers = {'o','v','diamond'};

resin_names = {{"3Dresyn Bioflex","A10 MB - IPA"},
               {"3Dresyn Bioflex","A10 MB - UNW2"},
               {"Asiga DentaGUIDE"},
               {"Asiga DentaGUM"},
               {"Formlabs Silicone","40A - IPA"},
               {"Formlabs Silicone","40A - IPA/BuOAc"},
               {"Liqcreate","Bio-Med Clear"},
               {"Phrozen AquaGray 8K"}};

%% Fitting Soft Resins with a Yeoh Model
min_frac = 0.05;

% this cell array will store all of the model fitting results
SOFT_MATERIALS_MODELS = cell(5,1);
resin_number = 0;

rng(0); % setting rng seed for repeatability

parameter_output_data = {"Resin","Sterilization","C1","C2","C3"};
metric_ouput_data = {"Resin","Sterilization","PlateID","SpecimenID","Stress 50% [MPa]","Stress 100% [MPa]","Ultimate Tensile Strength [MPa]","Strain at Break [%]"};

for i=[1,2,4,5,6] % loop through only the soft resins
    resin_number = resin_number + 1; % index into the SOFT_MATERIALS_MODEL cell array

    fprintf("\nResin: %s\n",resins{i})

    figure("Position",[100,100,1800,900],"Color","w"); hold all

    % storing the fits and statistics for the different sterilization
    % techniques
    fits = cell(3,1);
    stats = cell(3,1);
    

    for j=1:length(sterilization)
        col = colors{j}; %rand(3,1);

        % Compression Data
        inds = ismember(compress_data.Resin,resins(i)) & ismember(compress_data.Sterilization,sterilization(j));
        subset = compress_data(inds,:);
        min_frac = 0.05;
        
        % loop through all three plates and both samples per plate

        LAM_tot = []; % fitting all of the data simultaneously
        SIG_tot = [];
        indicator = [];

        YoungsMod = [];         % Young's modulus [MPa]
        strain_at_break = [];   % strain at failure [mm/mm]
        UTS = [];               % utilimate tensile strength [MPa]
        stress_50 = [];         % stress at 50% strain [MPa]
        stress_100 = [];        % stress at 100% strain [MPa]

        samp_count = 1;

        plates = unique(subset.PlateID);
        for k=1:length(plates)
            specimen = unique(subset{ismember(subset.PlateID,plates{k}),"SpecimenID"});
            for l=1:length(specimen)

                samp_ind = ismember(subset.PlateID,plates{k,1}) & ismember(subset.SpecimenID,specimen{l,1});
                samp_data = subset(samp_ind,:);

                % getting the rest geometry of the current sample
                geom_ind = ismember(geom_data.TestType,"Compression") & ...
                           ismember(geom_data.Resin,resins(i)) & ...
                           ismember(geom_data.Sterilization,sterilization(j)) & ...
                           ismember(geom_data.PlateID,plates{k}) & ...
                           ismember(geom_data.SpecimenID,specimen{l});

                use_sample = geom_data{geom_ind,"UseSample_NoForIfExperimentalError_"};
                sample_broke = geom_data{geom_ind,"SampleBroke_0Or1_"};

                if use_sample
                    L0 = geom_data{geom_ind,"ThicknessOrGageLength_mm_"};
                    A0 = geom_data{geom_ind,"Area_mm2_"};
                    
                    % pulling the data for the current sample
                    t = samp_data.Time_sec_;
                    dL = samp_data.Crosshead_mm_; 
                    F = samp_data.Load_N_;
                    
                    % plotting the raw force displacement data
                    subplot(length(sterilization),4,1 + 4*(j-1)); hold all
                    plot(-dL,-F,'.',"Color",col)
                    
                    % calculating the engineering stress
                    sig = F/A0;
        
                    % removing portions where slip occurs in rigid samples
                    ind_to_consider = (sig > min_frac*max(sig)); % only get stress above 5% of max to avoid region with slip in rigid
                    dL = dL(ind_to_consider); dL = dL - dL(1);
                    sig = sig(ind_to_consider);
                    % removing portions after peak force (plastic region)
                    [~,peak] = max(sig);
                    dL(peak:end) = [];
                    sig(peak:end) = [];
                    
                    sig = -1*sig; % negative because of compression
                    
                    % calculating the nominal stretch
                    lam = 1 - dL/L0;
                    
                    % using first order approximation to find the offset in the
                    % displacement
                    p1 = polyfit(lam(1:5),sig(1:5),1);
                    delL = L0*(p1(2)/p1(1) + 1); % offset found from linear fit
    
                    % interpolating the data to ensure equal weight between
                    % tensile and compression tests
                    lam_ = LAM(-dL,delL,L0);
                    lam_interp = linspace(min(lam_),max(lam_),N_interp_points);
                    sig_interp = interp1(lam_,sig,lam_interp);
                    
                    % plotting the current test
                    subplot(length(sterilization),4,2 + 4*(j-1)); hold all
                    plot(lam_,sig,'.',"Color",col)
                    
                    % storing the tensile data for fitting
                    LAM_tot = [LAM_tot, lam_interp];
                    SIG_tot = [SIG_tot, sig_interp];
                    indicator = [indicator, 0*lam_interp + samp_count];
    
                    samp_count = samp_count + 1;
                end
            

            end
        end

        % Tensile Data
        inds = ismember(tensile_data.Resin,resins(i)) & ismember(tensile_data.Sterilization,sterilization(j));
        subset = tensile_data(inds,:);
        plates = unique(subset.PlateID);
        for k=1:length(plates)
            
            specimen = unique(subset{ismember(subset.PlateID,plates{k}),"SpecimenID"});

            for l=1:length(specimen)
                
                % getting the data for the current sample
                samp_ind = ismember(subset.PlateID,plates{k,1}) & ismember(subset.SpecimenID,specimen{l,1});                
                samp_data = subset(samp_ind,:);
                
                % getting the rest geometry of the current sample
                geom_ind = ismember(geom_data.TestType,"Tensile") & ...
                           ismember(geom_data.Resin,resins(i)) & ...
                           ismember(geom_data.Sterilization,sterilization(j)) & ...
                           ismember(geom_data.PlateID,plates{k}) & ...
                           ismember(geom_data.SpecimenID,specimen{l});

                use_sample = geom_data{geom_ind,"UseSample_NoForIfExperimentalError_"};
                sample_broke = geom_data{geom_ind,"SampleBroke_0Or1_"};

                if use_sample
                    L0 = geom_data{geom_ind,"ThicknessOrGageLength_mm_"};
                    A0 = geom_data{geom_ind,"Area_mm2_"};
    
                    % pulling the data for the current sample
                    t = samp_data.Time_sec_;
                    dL = samp_data.Crosshead_mm_; 
                    F = samp_data.Load_N_;
                    
                    % plotting the raw force displacement data
                    subplot(length(sterilization),4,1 + 4*(j-1)); hold all
                    plot(dL,F,'.',"Color",col)
                    
                    % calculating the engineering stress
                    sig = F/A0;
        
                    % removing portions where slip occurs in rigid samples
                    ind_to_consider = (sig > min_frac*max(sig)); % only get stress above 5% of max to avoid region with slip in rigid
                    dL = dL(ind_to_consider); dL = dL - dL(1);
                    sig = sig(ind_to_consider);
                    % removing portions after peak force (plastic region)
                    [~,peak] = max(sig);
                    dL(peak:end) = [];
                    sig(peak:end) = [];
                    
                    % calculating the nominal stretch
                    lam = 1 + dL/L0;
                    
                    % using first order approximation to find the offset in the
                    % displacement
                    p1 = polyfit(lam(1:5),sig(1:5),1);
                    delL = L0*(p1(2)/p1(1) + 1); % offset found from linear fit
    
                    % interpolating the data to ensure equal weight between
                    % tensile and compression tests
                    lam_ = LAM(dL,delL,L0);
                    lam_interp = linspace(min(lam_),max(lam_),N_interp_points);
                    sig_interp = interp1(lam_,sig,lam_interp);
                    
                    % plotting the current test
                    subplot(length(sterilization),4,2 + 4*(j-1)); hold all
                    plot(lam_,sig,'.',"Color",col)
                    
                    if sample_broke
                        UTS_ = max(sig_interp);
                        elongation_at_break_ = 100*(max(lam_interp) - 1);
                    else
                        UTS_ = NaN;
                        elongation_at_break_ = NaN;
                    end

                    % storing the tensile data for fitting
                    LAM_tot = [LAM_tot, lam_interp];
                    SIG_tot = [SIG_tot, sig_interp];
                    indicator = [indicator, 0*lam_interp - samp_count];
                    
                    stress_50 = [stress_50, interp1(lam_,sig,1.5)];
                    stress_100 = [stress_100, interp1(lam_,sig,2)];
                    
                    UTS = [UTS; UTS_];
                    strain_at_break = [strain_at_break; elongation_at_break_];

                    metric_ouput_data = [metric_ouput_data; {resins{i},sterilization{j},plates{k},specimen{l},interp1(lam_,sig,1.5),interp1(lam_,sig,2), UTS_, elongation_at_break_}];
    
                    samp_count = samp_count + 1;
                end
            end
        end
   
        % plotting all of the data together
        sgtitle(resins(i))
        subplot(length(sterilization),4,3 + 4*(j-1)); hold all
        plot(LAM_tot,SIG_tot,'.','Color',col,'MarkerSize',8)

        % fitting the data to a Yeoh model
        fprintf("Sterilization: %s\n",sterilization{j})
        [fit_i,stat_i] = Bootstrap_Model_Fit(LAM_tot,SIG_tot,indicator,0);

        [C1_mean,C1_std] = RoundToSigFigs(mean(fit_i.C1),std(fit_i.C1));
        [C2_mean,C2_std] = RoundToSigFigs(mean(fit_i.C2),std(fit_i.C2));
        [C3_mean,C3_std] = RoundToSigFigs(mean(fit_i.C3),std(fit_i.C3));

        fprintf("\tYeoh Model Parameters: \n\t\tC1 = %s +/- %s MPa,\n\t\tC2 = %s +/- %s MPa, \n\t\tC3 = %s +/- %s MPa\n",num2str(C1_mean),num2str(C1_std),num2str(C2_mean),num2str(C2_std),num2str(C3_mean),num2str(C3_std))

        fit_i.strain_at_break = strain_at_break;
        fit_i.UTS = UTS;
        fit_i.E = 6*fit_i.C1;
%         stat_i.strain_at_break_mean = mean(100*strain_at_break);
%         stat_i.strain_at_break_std = std(100*strain_at_break);
%         stat_i.UTS_mean = mean(UTS);
%         stat_i.UTS_std = std(UTS);
        fit_i.stress_50 = stress_50;
        fit_i.stress_100 = stress_100;
        stat_i.stress_50_mean = mean(stress_50);
        stat_i.stress_50_std = std(stress_50);
        stat_i.stress_100_mean = mean(stress_100);
        stat_i.stress_100_std = std(stress_100);

        [E_mean,E_std] = RoundToSigFigs(6*mean(fit_i.C1),6*std(fit_i.C1));
        [s50_mean,s50_std] = RoundToSigFigs(mean(stress_50),std(stress_50));
        [s100_mean,s100_std] = RoundToSigFigs(mean(stress_100),std(stress_100));
        

        fprintf("\tStress-Strain Measures: \n\t\tYoung's Modulus: %s +/- %s MPa,\n\t\tStress at 0.5 Strain: %s +/- %s MPa,\n\t\tStress at 1.0 Strain: %s +/- %s MPa\n",num2str(E_mean),num2str(E_std),num2str(s50_mean),num2str(s50_std),num2str(s100_mean),num2str(s100_std))

        fits{j} = fit_i;
        stats{j} = stat_i;

        cur_out = cell(length(fit_i.C1),5);
        cur_out(:,1) = resins(i);
        cur_out(:,2) = sterilization(j);
        cur_out(:,3) = num2cell(fit_i.C1);
        cur_out(:,4) = num2cell(fit_i.C2);
        cur_out(:,5) = num2cell(fit_i.C3);

        parameter_output_data = [parameter_output_data; cur_out];
        
        
        % plotting the resulting model fits
        fill([stat_i.lam_all,fliplr(stat_i.lam_all)],[stat_i.sig_mean + stat_i.sig_std, fliplr(stat_i.sig_mean - stat_i.sig_std)],col*0.6,'FaceAlpha',0.2,'EdgeColor','none',"HandleVisibility","off");
        plot(stat_i.lam_all,stat_i.sig_mean,'--','color',col*0.2,"DisplayName","Yeoh Model")
        
    
        % annotating the figures
        legend("Location","NW")
        xlabel('Stretch [mm/mm]')
        ylabel('Engineering Stress [MPa]')
        subplot(length(sterilization),4,1 + 4*(j-1));
        xlabel("Displacement [mm]")
        ylabel("Force [N]")
    
        subplot(length(sterilization),4,2 + 4*(j-1));
        xlabel("Stretch [mm/mm]")
        ylabel("Engineering Stress [MPa]")
        title(sterilization(j))
    
    end

    % plotting the mean and standard deviation fits for each sterilization
    % technique
    subplot(length(sterilization),4,[4:4:12]); hold all
    for k=1:length(sterilization)
        fill([stats{k}.lam_all,fliplr(stats{k}.lam_all)],[stats{k}.sig_mean + stats{k}.sig_std, fliplr(stats{k}.sig_mean - stats{k}.sig_std)],colors{k}*0.6,'FaceAlpha',0.2,'EdgeColor','none');
        plot(stats{k}.lam_all,stats{k}.sig_mean,'--','color',colors{k}*0.2)
    end
    legend("Location","NW")
    xlabel("Stretch [mm/mm]")
    ylabel("Engineering Stress [MPa]")

    % Storing the fits and statistics for the current resin
    SOFT_MATERIALS_MODELS{resin_number} = { resins{i}, fits, stats };

end

writecell(parameter_output_data,out_file,"Sheet","Soft Resin Model Parameters");
writecell(metric_ouput_data,out_file,"Sheet","Soft Resin Stress Metrics");
save(saveDir+"SoftMaterialModels.mat","SOFT_MATERIALS_MODELS");

%% Soft Resin Summary Figures

ylims = {[-0.5,0.9],[-0.5,0.5],[-1,5],[-1,3.5],[-1,3.5]};
xlims = {[0.5,2.7],[0.5,2],[0.7,2.5],[0.6,3],[0.6,3]};

softresin_numbers = [1,2,4,5,6];

for i=1:5
    
    fits = SOFT_MATERIALS_MODELS{i}{2};
    stats = SOFT_MATERIALS_MODELS{i}{3};

    figure("Position",[100,100,1400,900],"Color","w");
    
    for k=1:length(sterilization)
        subplot(5,3,[k,k+3]); hold all

        indicators = fits{k}.indicator;
        indicator = unique(indicators);

        for j=1:length(indicator)
            ind = indicators == indicator(j);
            plot(fits{k}.lam(ind),fits{k}.sig(ind),'.',"Color",colors{k},"DisplayName","Data")
        end
        plot(stats{k}.lam_all,stats{k}.sig_mean,lines{k},"Color",'k',"DisplayName","Yeoh Model","LineWidth",1)
        xlabel("Stretch [mm/mm]")
        ylabel("Engineering Stress [MPa]")
        text(1,max(stats{k}.sig_mean),"R^2 = "+num2str(stats{k}.R2,3),"FontSize",12)
        ylim(ylims{i})
        xlim(xlims{i})
        
        set(gca,"FontSize",15)

    end
        set(gca,"FontSize",15)

    
    subplot(5,3,6+[1,2,4,5,7,8]); hold all
    for k=1:length(sterilization)
        fill([stats{k}.lam_all,fliplr(stats{k}.lam_all)],[stats{k}.sig_mean + stats{k}.sig_std, fliplr(stats{k}.sig_mean - stats{k}.sig_std)],colors{k},'FaceAlpha',alpha,'EdgeColor','none',"HandleVisibility","off");
        plot(stats{k}.lam_all,stats{k}.sig_mean,lines{k},'color',colors{k},"DisplayName",sterilization{k},"LineWidth",1.5)
    end
    ylim(ylims{i})
    xlim(xlims{i})
    xlabel("Stretch [mm/mm]")
    ylabel("Engineering Stress [MPa]")
    legend("Location","northwest")
%     xlim([0.5,3])
%     ylim([-1,5])
    set(gca,"FontSize",15)


    subplot(5,3,9); hold all
    for k=1:length(sterilization)
        bar(k,stats{k}.C1_mean,0.5,"FaceColor",colors{k})
        errorbar(k,stats{k}.C1_mean,stats{k}.C1_std,"vertical","Marker","none","Color","k","MarkerFaceColor",colors{k},"linewidth",1.5)
    end
    
%     xlabel("Sterilization")
    ylabel("C_1 [MPa]")
    xlim([0.5,3.5])
    
    xticks([])
%     xticklabels(sterilization)
    set(gca,"FontSize",15)


    subplot(5,3,12); hold all
    for k=1:length(sterilization)
        bar(k,stats{k}.C2_mean,0.5,"FaceColor",colors{k})
        errorbar(k,stats{k}.C2_mean,stats{k}.C2_std,"vertical","Marker","none","Color","k","MarkerFaceColor",colors{k},"linewidth",1.5)
    end
    
%     xlabel("Sterilization")
    ylabel("C_2 [MPa]")
    xlim([0.5,3.5])
    xticks([])
%     xticklabels(sterilization)
    set(gca,"FontSize",15)

    subplot(5,3,15); hold all
    for k=1:length(sterilization)
        bar(k,stats{k}.C3_mean,0.5,"FaceColor",colors{k})
        errorbar(k,stats{k}.C3_mean,stats{k}.C3_std,"vertical","Marker","none","Color","k","MarkerFaceColor",colors{k},"linewidth",1.5)
    end
    
    xlabel("Sterilization")
    ylabel("C_3 [MPa]")
    xlim([0.5,3.5])
    xticks(1:3)
    xticklabels(sterilization)
    sgtitle(resin_names{softresin_numbers(i)},"FontSize",15); %upper(SOFT_MATERIALS_MODELS{i}{1}))
    set(gca,"FontSize",15)
    
    saveas(gcf,saveDir+"Output Figures\"+SOFT_MATERIALS_MODELS{i}{1}+" Model Summary.svg")


    figure("Position",[100,100,250,700],"Color","w"); hold all
    subplot(3,1,1); hold all
    for m=3:-1:1
        Es = stats{m}.C1_mean*6;
        dEs = stats{m}.C1_std*6;
        bar(m,mean(Es),0.8,"FaceColor",colors{m})
        errorbar(m,Es,dEs,"vertical","-k","MarkerFaceColor",'k',"LineWidth",1.5);

        if m==1
            plot([0.5,3.5],mean(Es)*[1.05,1.05],'--k')
            plot([0.5,3.5],mean(Es)*[0.95,0.95],'--k')
        end
    end
%     xlabel("Sterilization")
    ylabel({"Young's","Modulus [MPa]"})
    xlim([0.5,3.5])
    ylim([0,2.5]*(0.5 + 0.5*(SOFT_MATERIALS_MODELS{i}{1} == "asiga dentagum")))
%     xticks(1:3)
%     xticklabels(sterilization)
    xticks([])
    set(gca,"FontSize",12)


    subplot(3,1,2); hold all
    for m=3:-1:1
        stress_50 = fits{m}.stress_50;
        bar(m,mean(stress_50),0.8,"FaceColor",colors{m})
        errorbar(m,mean(stress_50),std(stress_50),"vertical","-k","MarkerFaceColor",'k',"LineWidth",1.5);
        plot(m - 0.25 + 0*stress_50,stress_50,markers{m},"Color","k","MarkerFaceColor",colors{m},"MarkerSize",3);

        if m==1
            plot([0.5,3.5],mean(stress_50)*[1.05,1.05],'--k')
            plot([0.5,3.5],mean(stress_50)*[0.95,0.95],'--k')
        end
    end
%     xlabel("Sterilization")
    ylabel({"Stress at 50%","Strain [MPa]"})
    xlim([0.5,3.5])
    ylim([0,1.1]*(0.4 + 0.6*(SOFT_MATERIALS_MODELS{i}{1} == "asiga dentagum")))
%     xticks(1:3)
%     xticklabels(sterilization)
    xticks([])
    set(gca,"FontSize",12)

    subplot(3,1,3); hold all
    for m=3:-1:1
        stress_100 = fits{m}.stress_100;
        bar(m,mean(stress_100),0.8,"FaceColor",colors{m})
        errorbar(m,mean(stress_100),std(stress_100),"vertical","-k","MarkerFaceColor",'k',"LineWidth",1.5);
        plot(m - 0.25 + 0*stress_100,stress_100,markers{m},"Color","k","MarkerFaceColor",colors{m},"MarkerSize",3);

        if m==1
            plot([0.5,3.5],mean(stress_100)*[1.05,1.05],'--k')
            plot([0.5,3.5],mean(stress_100)*[0.95,0.95],'--k')
        end
    end
    xlabel("Sterilization")
    ylabel({"Stress at 100%","Strain [MPa]"})
    xlim([0.5,3.5])
    ylim([0,2.5]*(0.45 + 0.55*(SOFT_MATERIALS_MODELS{i}{1} == "asiga dentagum")))
    xticks(1:3)
    xticklabels(sterilization)
    set(gca,"FontSize",12)

    sgtitle(resin_names{softresin_numbers(i)},"FontSize",12); %upper(SOFT_MATERIALS_MODELS{i}{1}))

    saveas(gcf,saveDir+"Output Figures\"+SOFT_MATERIALS_MODELS{i}{1}+" Parameters.svg")

    
    if ~((SOFT_MATERIALS_MODELS{i}{1} == "formlabs silicone mix") | (SOFT_MATERIALS_MODELS{i}{1} == "formlabs silicone ipa"))
        figure("Position",[100,100,250,500],"Color","w"); hold all
        subplot(2,1,1); hold all
        for m=3:-1:1
            UTS = fits{m}.UTS;
            bar(m,mean(UTS),0.8,"FaceColor",colors{m})
            errorbar(m,mean(UTS),std(UTS),"vertical","-k","MarkerFaceColor",'k',"LineWidth",1.5);
            plot(m - 0.25 + 0*UTS,UTS,markers{m},"Color","k","MarkerFaceColor",colors{m},"MarkerSize",3);

            if m==1
                plot([0.5,3.5],mean(UTS)*[1.05,1.05],'--k')
                plot([0.5,3.5],mean(UTS)*[0.95,0.95],'--k')
            end
        end

        xlabel("Sterilization")
        ylabel({"Ultimate Tensile","Strength [MPa]"})
        xlim([0.5,3.5])
        ylim([0,0.7]*(1 + 5*(SOFT_MATERIALS_MODELS{i}{1} == "asiga dentagum")))
        xticks(1:3)
        xticklabels(sterilization)
        set(gca,"FontSize",12)


        subplot(2,1,2); hold all
        for m=3:-1:1
            SAB = fits{m}.strain_at_break;
            bar(m,mean(SAB),0.8,"FaceColor",colors{m})
            errorbar(m,mean(SAB),std(SAB),"vertical","-k","MarkerFaceColor",'k',"LineWidth",1.5);
            plot(m - 0.25 + 0*SAB,SAB,markers{m},"Color","k","MarkerFaceColor",colors{m},"MarkerSize",3);

            if m==1
                plot([0.5,3.5],mean(SAB)*[1.05,1.05],'--k')
                plot([0.5,3.5],mean(SAB)*[0.95,0.95],'--k')
            end
        end


        subplot(2,1,2); hold all
        xlabel("Sterilization")
        ylabel({"Elongation","at Break [%]"})
        xlim([0.5,3.5])
        ylim([0,170])
        xticks(1:3)
        xticklabels(sterilization)
        set(gca,"FontSize",12)

        sgtitle(resin_names{softresin_numbers(i)},"FontSize",12); %upper(SOFT_MATERIALS_MODELS{i}{1}))

        saveas(gcf,saveDir+"Output Figures\"+SOFT_MATERIALS_MODELS{i}{1}+" FailureStats.svg")
    end


end

%% Calculating percent differences between sterilization groups

metrics = {"Young's Modulus","Stress at 50%%","Stress at 100%%","Elongation at Break","UTS"};
names = {"E","stress_50","stress_100","strain_at_break","UTS"};
disp("Percent Change in the Mean:")
for i=1:5
    % for all resin groups
    data_i = SOFT_MATERIALS_MODELS{i};
    disp(data_i{1})
    for j=1:5
        % for all metrics
        fprintf("\t\t Metric: "+metrics{j}+"\n")
        
        % non sterile data
        ns_data = data_i{2}{1};
        ns_mean = mean(ns_data.(names{j}));
        ns_std = std(ns_data.(names{j}));

        % autoclave data
        ac_data = data_i{2}{2};
        ac_mean = mean(ac_data.(names{j}));
        ac_std = std(ac_data.(names{j}));

        % ethanol data
        et_data = data_i{2}{3};
        et_mean = mean(et_data.(names{j}));
        et_std = std(et_data.(names{j}));
        
        % Autoclave difference
        ac_diff = 100*(ac_mean - ns_mean)/ns_mean;
        ac_diff_err = abs(ac_diff)*sqrt( ( ac_std / ac_mean )^2 + ( ns_std / (ns_mean) )^2);
        [ac_diff,ac_diff_err] = RoundToSigFigs(ac_diff,ac_diff_err);

        fprintf("\t\t\tNonsterile value: %s +/- %s\n",num2str(ns_mean),num2str(ns_std))
        fprintf("\t\t\tAutoclave value: %s +/- %s\n",num2str(ac_mean),num2str(ac_std))        
        
        fprintf("\t\t\t\tAutoclave Difference: %s +/- %s %%\n",num2str(ac_diff),num2str(ac_diff_err))

        % Ethanol difference
        et_diff = 100*(et_mean - ns_mean)/ns_mean;
        et_diff_err = abs(et_diff)*sqrt( ( et_std / et_mean )^2 + ( ns_std / (ns_mean) )^2);
        [et_diff,et_diff_err] = RoundToSigFigs(et_diff,et_diff_err);

        fprintf("\t\t\tEthanol value: %s +/- %s\n",num2str(et_mean),num2str(et_std))        
        
        fprintf("\t\t\t\tEthanol Difference: %s +/- %s %%\n\n\n",num2str(et_diff),num2str(et_diff_err))

    
    end

end



%% Fitting Rigid Resins with a Linear Elasticity Model
min_frac = 0.;
% colors = {[0.2,0.8,0.6], [0.8,0.2,0.5], [0.2,0.6,0.8]};


peak_stress = zeros(8,1);
peak_stress(3) = 60;
peak_stress(7) = 60;
peak_stress(8) = 30;

peak_c_stress = zeros(8,1);
peak_c_stress(3) = 11;
peak_c_stress(7) = 11;
peak_c_stress(8) = 14;


% this cell array will store all of the model fitting results
RIGID_MATERIALS_MODELS = cell(3,1);
resin_number = 0;

rigid_tensile_out_data = {"Resin","Sterilization","PlateID","SpecimenID","Young's Modulus [MPa]","Ultimate Tensile Strength [MPa]","Elongation at Break [%]"};
rigid_compressive_out_data = {"Resin","Sterilization","PlateID","SpecimenID","Compressive Modulus [MPa]"};

for i=[3,7,8] % loop through only the rigid resins
    resin_number = resin_number + 1; % index into the SOFT_MATERIALS_MODEL cell array

    fprintf("\nResin: %s\n",resins{i})

    figure("Position",[100,100,900,520],"Color","w"); hold all

    % storing the fits and statistics for the different sterilization
    % techniques
    fits = cell(3,1);
    stats = cell(3,1);

    for j=1:length(sterilization)
        col = colors{j}; %rand(3,1);

        % Compression Data
        inds = ismember(compress_data.Resin,resins(i)) & ismember(compress_data.Sterilization,sterilization(j));
        subset = compress_data(inds,:);
        
        % loop through all three plates and both samples per plate
        YoungsMod = [];         % Young's modulus [MPa]
        CompMod = [];           % Compressive Modulus [MPa]
        strain_at_break = [];   % strain at failure [mm/mm]
        UTS = [];               % utilimate tensile strength [MPa]


        samp_count = 1;
        
        subplot(2,3,j+3); hold all

        plates = unique(subset.PlateID);
        for k=1:length(plates)
            specimen = unique(subset{ismember(subset.PlateID,plates{k}),"SpecimenID"});
            for l=1:length(specimen)

                samp_ind = ismember(subset.PlateID,plates{k,1}) & ismember(subset.SpecimenID,specimen{l,1});
                samp_data = subset(samp_ind,:);

                % getting the rest geometry of the current sample
                geom_ind = ismember(geom_data.TestType,"Compression") & ...
                           ismember(geom_data.Resin,resins(i)) & ...
                           ismember(geom_data.Sterilization,sterilization(j)) & ...
                           ismember(geom_data.PlateID,plates{k}) & ...
                           ismember(geom_data.SpecimenID,specimen{l});

                use_sample = geom_data{geom_ind,"UseSample_NoForIfExperimentalError_"};
                sample_broke = geom_data{geom_ind,"SampleBroke_0Or1_"};

                if use_sample
                    L0 = geom_data{geom_ind,"ThicknessOrGageLength_mm_"};
                    A0 = geom_data{geom_ind,"Area_mm2_"};
                    
                    % pulling the data for the current sample
                    t = samp_data.Time_sec_;
                    dL = samp_data.Crosshead_mm_; 
                    F = samp_data.Load_N_;
                    
                    % displacement rate
                    if (i==3)&(j==1)&(k==1)&(l==1)
                        % only show for first sample
                        dLdt = 100*mean(diff(dL/L0)./diff(t))*60;
                        fprintf("Test Rate for Compression: %.3f %%/min\n",dLdt)
                    end

                    % calculating the engineering stress
                    sig = F/A0;
                    dL = dL - dL(1);

                    sig = 1*sig; % negative because of compression
                    
                    % calculating the nominal stretch
                    lam = 1 + dL/L0;
                    
                    % interpolating the data to ensure equal weight between
                    % tensile and compression tests
                    lam_ = LAM(dL,0,L0);
                    lam_interp = linspace(min(lam_),max(lam_),N_interp_points);
                    sig_interp = interp1(lam_,sig,lam_interp);
    
                    e1 = diff(sig_interp)./diff(lam_interp);
                    e1 = smooth(e1,10,"moving");
                    
                    [~,ind_] = max(abs(e1));
    
                    if ind_+1 > length(lam_interp)
                        linear_region = length(lam_interp)-2:length(lam_interp);
                    elseif ind_<2
                        linear_region = 1:3;
                    else
                        linear_region = ind_-1:ind_+1;
                    end
    
                    eps_linear = lam_interp(linear_region)' - 1;
                    sig_linear = sig_interp(linear_region)';
    
                    p1 = polyfit(eps_linear,sig_linear,1);
    
                    E1 = p1(1);
                    delL = L0*(p1(2)/p1(1));
                    lam_ = LAM(dL,delL,L0);
                    lam_interp = linspace(min(lam_),max(lam_),N_interp_points);
                    
                    % plotting the current test
                    plot(100*(lam_interp-1),sig_interp,".","Color",col)
                    plot(100*([1,max(lam_)]-1),E1*([1,max(lam_)]-1),'--',"Color","k")

                    xlabel("Compressive Strain [%]")
                    ylabel({"Compressive Engineering","Stress [MPa]"})
                    set(gca,"FontSize",12)

                    xlim(100*([0.97,1.28]-1))
                    ylim([0,peak_c_stress(i)])

                    CompMod = [CompMod, E1];

                    rigid_compressive_out_data = [rigid_compressive_out_data; {resins{i},sterilization{j},plates{k},specimen{l},E1}];
    
                    samp_count = samp_count + 1;
                end

            end
        end

        subplot(2,3,j); hold all
        % Tensile Data
        inds = ismember(tensile_data.Resin,resins(i)) & ismember(tensile_data.Sterilization,sterilization(j));
        subset = tensile_data(inds,:);
        plates = unique(subset.PlateID);
        for k=1:length(plates)
            
            specimen = unique(subset{ismember(subset.PlateID,plates{k}),"SpecimenID"});

            for l=1:length(specimen)
                
                % getting the data for the current sample
                samp_ind = ismember(subset.PlateID,plates{k,1}) & ismember(subset.SpecimenID,specimen{l,1});                
                samp_data = subset(samp_ind,:);
                
                % getting the rest geometry of the current sample
                geom_ind = ismember(geom_data.TestType,"Tensile") & ...
                           ismember(geom_data.Resin,resins(i)) & ...
                           ismember(geom_data.Sterilization,sterilization(j)) & ...
                           ismember(geom_data.PlateID,plates{k}) & ...
                           ismember(geom_data.SpecimenID,specimen{l});

                use_sample = geom_data{geom_ind,"UseSample_NoForIfExperimentalError_"};
                sample_broke = geom_data{geom_ind,"SampleBroke_0Or1_"};

                if use_sample
                    L0 = geom_data{geom_ind,"ThicknessOrGageLength_mm_"};
                    A0 = geom_data{geom_ind,"Area_mm2_"};
    
                    % pulling the data for the current sample
                    t = samp_data.Time_sec_;
                    dL = samp_data.Crosshead_mm_; 
                    F = samp_data.Load_N_;

                     % displacement rate
                    if (i==7)&(j==1)&(k==1)&(l==1)
                        % only show for first sample
                        dLdt = 100*mean(diff(dL/L0)./diff(t))*60;
                        fprintf("Test Rate for Tension: %.3f %%/min\n",dLdt)
                    end

                    % calculating the engineering stress
                    sig = F/A0;
                    dL = dL - dL(1);
                    
                    % calculating the nominal stretch
                    lam = 1 + dL/L0;
                    
                    % interpolating the data to ensure equal weight between
                    % tensile and compression tests
                    lam_ = LAM(dL,0,L0);
                    lam_interp = linspace(min(lam_),max(lam_),50);
                    sig_interp = interp1(lam_,sig,lam_interp);
                    
                    e1 = diff(sig_interp)./diff(lam_interp);
                    e1 = smooth(e1,5,"moving");
                    
                    [~,ind_] = max(e1);
                    if ind_ < 2
                        % highest slope is right at the beginning
                        linear_region = 1:3;
                    else
                        linear_region = ind_-1:ind_+1;
                    end
                    eps_linear = lam_interp(linear_region)' - 1;
                    sig_linear = sig_interp(linear_region)';
    
                    p1 = polyfit(eps_linear,sig_linear,1);
    
                    E1 = p1(1);
                    delL = L0*(p1(2)/p1(1));
                    lam_ = LAM(dL,delL,L0);
                    lam_interp = linspace(min(lam_),max(lam_),50);
    
    
                    % plotting the current test
                    plot(100*(lam_interp-1),sig_interp,".","Color",col)
                    plot(100*([1,max(lam_)]-1),E1*([1,max(lam_)]-1),'--',"Color","k")
    
                   samp_count = samp_count + 1;
                    
                   YoungsMod = [YoungsMod, E1];
                   strain_at_break = [strain_at_break, max(lam_interp)-1];
                   UTS = [UTS, max(sig_interp)];

                   rigid_tensile_out_data = [rigid_tensile_out_data; {resins{i},sterilization{j},plates{k},specimen{l},E1,max(sig_interp),100*(max(lam_interp)-1)}];

                end
            end
        end


               
        fits{j} = {YoungsMod,100*strain_at_break,UTS,CompMod};
        % annotating the figures
        xlabel('Tensile Strain [%]')
        ylabel({"Tensile Engineering","Stress [MPa]"})
        xlim(100*([0.95,1.4]-1))
        ylim([0,peak_stress(i)])
        sgtitle(resin_names{i},"FontSize",12)

        set(gca,"FontSize",12)

        [E_mean,E_std] = RoundToSigFigs(mean(YoungsMod),std(YoungsMod));
        [C_mean,C_std] = RoundToSigFigs(mean(CompMod),std(CompMod));
        [SAB_mean,SAB_std] = RoundToSigFigs(mean(100*strain_at_break),std(100*strain_at_break));
        [UTS_mean,UTS_std] = RoundToSigFigs(mean(UTS),std(UTS));
        
        

        fprintf("Sterilization: %s\n",sterilization{j})
        fprintf("\tYoung's Modulus: %s +/- %s MPa\n",num2str(E_mean),num2str(E_std));
        fprintf("\tCompressive Modulus: %s +/- %s MPa\n",num2str(C_mean),num2str(C_std));
        fprintf("\tStrain at Break: %s +/- %s [percent]\n",num2str(SAB_mean),num2str(SAB_std));
        fprintf("\tUltimate Tensile Strength: %s +/- %s MPa\n",num2str(UTS_mean),num2str(UTS_std));
        
        

    
    end
    
    RIGID_MATERIALS_MODELS{resin_number} = [resins(i);fits];

    saveas(gcf,saveDir+"Output Figures\"+RIGID_MATERIALS_MODELS{resin_number}{1}+" StressStrainData.svg")


    figure("Position",[100,100,250,1200],"Color","w"); hold all
    subplot(4,1,1); hold all
    for m=3:-1:1
        Es = fits{m}{1};
        bar(m,mean(Es),0.8,"FaceColor",colors{m})
        errorbar(m,mean(Es),std(Es),"vertical","-k","MarkerFaceColor",'k',"LineWidth",1.5);
        plot(m - 0.25 + 0*Es,Es,markers{m},"Color","k","MarkerFaceColor",colors{m},"MarkerSize",3);

        if m==1
            plot([0.5,3.5],mean(Es)*[1.05,1.05],'--k')
            plot([0.5,3.5],mean(Es)*[0.95,0.95],'--k')
        end
    end
    xlabel("Sterilization")
    ylabel({"Young's","Modulus [MPa]"})
    xlim([0.5,3.5])
    xticks(1:3)
    ylim([0,850])
    xticklabels(sterilization)
    set(gca,"FontSize",12)


    subplot(4,1,3); hold all
    for m=3:-1:1
        UTSs = fits{m}{3};
        bar(m,mean(UTSs),0.8,"FaceColor",colors{m})
        errorbar(m,mean(UTSs),std(UTSs),"vertical","-k","MarkerFaceColor",'k',"LineWidth",1.5);
        plot(m - 0.25 + 0*UTSs,UTSs,markers{m},"Color","k","MarkerFaceColor",colors{m},"MarkerSize",3);

        if m==1
            plot([0.5,3.5],mean(UTSs)*[1.05,1.05],'--k')
            plot([0.5,3.5],mean(UTSs)*[0.95,0.95],'--k')
        end
    end
    xlabel("Sterilization")
    ylabel({"Utilimate Tensile","Strength [MPa]"})
    xlim([0.5,3.5])
    xticks(1:3)
    ylim([0,60])
    xticklabels(sterilization)
    set(gca,"FontSize",12)

    subplot(4,1,4); hold all
    for m=3:-1:1
        epss = fits{m}{2};
        bar(m,mean(epss),0.8,"FaceColor",colors{m})
        errorbar(m,mean(epss),std(epss),"vertical","-k","MarkerFaceColor",'k',"LineWidth",1.5);
        plot(m - 0.25 + 0*epss,epss,markers{m},"Color","k","MarkerFaceColor",colors{m},"MarkerSize",3);

        if m==1
            plot([0.5,3.5],mean(epss)*[1.05,1.05],'--k')
            plot([0.5,3.5],mean(epss)*[0.95,0.95],'--k')
        end
    end
    xlabel("Sterilization")
    ylabel({"Elongation","at Break [%]"})
    xlim([0.5,3.5])
    xticks(1:3)
    ylim([0,45])
    xticklabels(sterilization)
    set(gca,"FontSize",12)

    subplot(4,1,2); hold all
    for m=3:-1:1
        CM = fits{m}{4};
        bar(m,mean(CM),0.8,"FaceColor",colors{m})
        errorbar(m,mean(CM),std(CM),"vertical","-k","MarkerFaceColor",'k',"LineWidth",1.5);
        plot(m - 0.25 + 0*CM,CM,markers{m},"Color","k","MarkerFaceColor",colors{m},"MarkerSize",3);

        if m==1
            plot([0.5,3.5],mean(CM)*[1.05,1.05],'--k')
            plot([0.5,3.5],mean(CM)*[0.95,0.95],'--k')
        end
    end
    xlabel("Sterilization")
    ylabel({"Compressive"," Modulus [MPa]"})
    xlim([0.5,3.5])
    xticks(1:3)
    ylim([0,250])
    xticklabels(sterilization)

    sgtitle(resin_names{i},"FontSize",12)

    set(gca,"FontSize",12)

    saveas(gcf,saveDir+"Output Figures\"+RIGID_MATERIALS_MODELS{resin_number}{1}+" Properties.svg")



end

writecell(rigid_compressive_out_data,out_file,"Sheet","Rigid Resin Compressive Metrics");
writecell(rigid_tensile_out_data,out_file,"Sheet","Rigid Resin Tensile Metrics");
save(saveDir+"RigidMaterialModels.mat","RIGID_MATERIALS_MODELS");

%% calculating percent differences in the means of the different groups

metrics = {"Young's Modulus","Elongation at Break","UTS","Compressive Modulus"};
disp("Percent Change in the Mean:")
for i=1:3
    % for all resin groups
    data_i = RIGID_MATERIALS_MODELS{i};
    disp(data_i{1})
    for j=1:4
        % for all metrics
        fprintf("\t\t Metric: "+metrics{j}+"\n")
        
        % non sterile data
        ns_data = data_i{2};
        ns_mean = mean(ns_data{j});
        ns_std = std(ns_data{j});

        % autoclave data
        ac_data = data_i{3};
        ac_mean = mean(ac_data{j});
        ac_std = std(ac_data{j});

        % ethanol data
        et_data = data_i{4};
        et_mean = mean(et_data{j});
        et_std = std(et_data{j});
        
        % Autoclave difference
        ac_diff = 100*(ac_mean - ns_mean)/ns_mean;
        ac_diff_err = abs(ac_diff)*sqrt( ( ac_std / ac_mean )^2 + ( ns_std / (ns_mean) )^2);
        [ac_diff,ac_diff_err] = RoundToSigFigs(ac_diff,ac_diff_err);

        fprintf("\t\t\tNonsterile value: %s +/- %s\n",num2str(ns_mean),num2str(ns_std))
        fprintf("\t\t\tAutoclave value: %s +/- %s\n",num2str(ac_mean),num2str(ac_std))        
        
        fprintf("\t\t\t\tAutoclave Difference: %s +/- %s %%\n",num2str(ac_diff),num2str(ac_diff_err))

        % Ethanol difference
        et_diff = 100*(et_mean - ns_mean)/ns_mean;
        et_diff_err = abs(et_diff)*sqrt( ( et_std / et_mean )^2 + ( ns_std / (ns_mean) )^2);
        [et_diff,et_diff_err] = RoundToSigFigs(et_diff,et_diff_err);

        fprintf("\t\t\tEthanol value: %s +/- %s\n",num2str(et_mean),num2str(et_std))        
        
        fprintf("\t\t\t\tEthanol Difference: %s +/- %s %%\n\n\n",num2str(et_diff),num2str(et_diff_err))

    
    end

end

