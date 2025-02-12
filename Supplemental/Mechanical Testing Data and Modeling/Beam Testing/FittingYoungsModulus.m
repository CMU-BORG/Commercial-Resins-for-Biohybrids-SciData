%% Fitting Young's Modulus to Beam Data

clear all; close all; clc

addpath("..\")
to_save = 0;

resins = ["Phrozen","Liqcreate","DentaGuide"];
samp_inds = [1,2,3,4;2,3,4,5;0,0,0,0];

Youngs_paper = [384,709,676];
err_Youngs_paper = [37;50;60];

Youngs_company = [2640.5,2000,0];
err_Youngs_company = [207.5,0,0];

YoungsModuli = 0*samp_inds;
err_E = 0*samp_inds;

for j=1:2
    resin = resins(j);
    fprintf("Resin: "+resin+"\n")
    
    for i=1:size(samp_inds,2)
        if i<=2
            state = "control";
        else
            state = "soaked";
        end

        fprintf("\tSample: "+num2str(samp_inds(j,i))+"\n")
        %% Choosing the Sample
        sample = num2str(samp_inds(j,i));
        folder = resin+" Beams\Sample "+sample+"\";
        file = resin+"Sample"+sample+"_data.csv";
        
        rotate_scale = 1;
        
        %% Reading in the Data
        
        % Composite Image
        image = imread(folder+"Overlay.png");
        [u,v,~] = size(image);
        
        % Digitize data
        data = readmatrix(folder+file,"NumHeaderLines",2);
        
        if rotate_scale
            scale_x = data(1:2,1);
            scale_y = data(1:2,2);
            scale = sqrt(diff(scale_x).^2 + diff(scale_y).^2) / ((130-70)/1000);
        else
            scale = data(1:2,1); scale = diff(scale) / ((130-70)/1000);
        end
        
        % displacement data
        rest_x = data(:,3) / scale; rest_x = [rest_x(end);rest_x(1:end-1)];
        rest_y = data(:,4) / scale; rest_y = [rest_y(end);rest_y(1:end-1)];
        load_100g_x = data(:,5) / scale; load_100g_x = [load_100g_x(end);load_100g_x(1:end-1)];
        load_100g_y = data(:,6) / scale; load_100g_y = [load_100g_y(end);load_100g_y(1:end-1)];
        load_200g_x = data(:,7) / scale; load_200g_x = [load_200g_x(end);load_200g_x(1:end-1)];
        load_200g_y = data(:,8) / scale; load_200g_y = [load_200g_y(end);load_200g_y(1:end-1)];
        
        % re-referencing data
        y_origin = rest_y(1);
        rest_y = y_origin - rest_y;
        load_100g_y = y_origin - load_100g_y;
        load_200g_y = y_origin - load_200g_y;
        
        x_origin = rest_x(1);
        rest_x = -x_origin + rest_x;
        load_100g_x = -x_origin + load_100g_x;
        load_200g_x = -x_origin + load_200g_x;
        
        %% Fitting Spline to Rest Geometry
        n = 1000; 
        
        % reinterpolating the data at a higher resolution
        x_ = linspace(rest_x(1),rest_x(end),n); 
        y_ = spline(rest_x,rest_y,x_);
        
        % calculating the arclength parameters
        dx = diff(x_);
        dy = diff(y_);
        ds = [0,sqrt(dx.^2 + dy.^2)];
        
        L_ = sum(ds);
        s = cumsum(ds) / L_; % getting normalize arclength coordinate of points
        
        theta = atan2(dy,dx); theta = [theta(1),theta];
        kappa_0 = diff(theta)./diff(s); kappa_0 = [kappa_0(1),kappa_0];
        
        % reinterpolating the rest position
        n = 25; % after convergence test, there is less than a 2% error is the Young's Modulus approximation
        convergence_error = 0.04; % worst case error is approximated as 4%
        S = linspace(0,1,n);
        x0 = interp1(s,x_,S,"linear","extrap")';
        y0 = interp1(s,y_,S,"linear","extrap")';
        theta0 = interp1(s,theta,S,"linear","extrap")';
        
        S_data = zeros(length(rest_x),1);
        for k=1:length(S_data)
            [~,ind] = min(sqrt( (x0 - rest_x(k)).^2 + (y0 - rest_y(k)).^2 ));
            S_data(k) = S(ind);
        end
        
        %% Setting up the Beam Geometry
        force_scaler = 1;

        % geometry parameters
        w1 = 10e-3;                    % [m] width of beam
        w2 = 3e-3;                     % [m] thickness of beam
        A = w1*w2;                     % [m^2] cross-sectional area
        I = (1/12)*(w1 * w2^3);        % [m^4] polar moment of inertia
        L = L_;                         % [m] length of the beam
        
        ds = L/n;
        dt = 1e-8;
        rho = 1;
        g = 9.81 / force_scaler;
        mass = 0*A*L*rho;
        
        % finite difference matrix
        [dds,~] = D1_matrix(n,n,L,L);
        dds = ds*dds(1:n,1:n);
        
        % de-dimensionalizing parameters
        Nc = (dt^2) / (rho*A*ds^2);
        Mc = (dt^2) / (rho*I*ds);
        
        E = 2200e6 / force_scaler;
        params = struct();
        params.A = A;
        params.I = I;
        params.L = L;
        params.dt = dt;
        params.ds = ds;
        params.dds = dds;
        params.Mc = Mc;
        params.Nc = Nc;
        params.n = n;
        params.w = w2;
        params.DELTA_t = 2^-0;
        
        mass = 0.200;
        fy = -mass*g;
        % [x,y] = HingeDisplacement(E,x0,y0,theta0,params,fy,S_data);
        %%
        err = @(E) ModelError(E,x0,y0,theta0,params,S_data,0.100,0.200,...
                              load_100g_x,load_100g_y,load_200g_x,load_200g_y);
        
        E = fminsearch(err,E);%,optimset("display","iter"));
        
        [x1,y1] = HingeDisplacement(E,x0,y0,theta0,params,-0.1*g,S_data);
        [x2,y2] = HingeDisplacement(E,x0,y0,theta0,params,-0.2*g,S_data);
        
        %% Plotting Results
        E_ = round(force_scaler*E/(1e9),2)*1000;
        E_err = E_*convergence_error;
        YoungsModuli(j,i) = E_;
        err_E(j,i) = E_err;
        fprintf("\tYoung's Modulus: %.0f +/- %.0f MPa\n",E_,E_err)
        
        
        figure("Color","w"); 
        imshow(image,"XData",[0,v/scale] - x_origin,"YData",y_origin - [0,u/scale])
        set(gca,'YDir','normal')
        
        hold all
        plot(rest_x,rest_y,'or','MarkerFaceColor','r',"DisplayName","0 g")
        plot(load_100g_x,load_100g_y,'sg','MarkerFaceColor','g',"DisplayName","100 g")
        plot(load_200g_x,load_200g_y,'^b','MarkerFaceColor','b',"DisplayName","200 g")
        
        plot(x0,y0,'--k',"LineWidth",2,"DisplayName","Model Prediction")
        plot(x0,y0,'--r',"LineWidth",2,"HandleVisibility","off")
        plot(x1,y1,'--g',"LineWidth",2,"HandleVisibility","off")
        plot(x2,y2,'--b',"LineWidth",2,"HandleVisibility","off")
        
        legend("Location","NW")
        title(resin+" Sample "+sample + " | E = "+num2str(E_)+" MPa")
        set(gca,"FontSize",15)
        if to_save
            saveas(gcf,folder+resin+"_Sample_"+sample+".svg","svg")
        end
    
    end
end

%%
figure("Position",[100,100,1500,400],"color","w"); hold all

w_space = 0.15;
mark_size = 4;
line_size = 1;


for j=1:length(resins)
    subplot(1,length(resins),j); hold all
    
    % control
    mean_ = mean(YoungsModuli(j,1:2));
    std_ = std(YoungsModuli(j,1:2));
    err_ = (1/2)*sqrt(sum(err_E(j,1:2).^2));

    col = 0.5 + 0.5*rand(3,1);
    for i=1:2 %size(samp_inds,2)
        errorbar(2 - w_space + 2*w_space*((i-1)/(2-1)), YoungsModuli(j,i), err_E(j,i),"vertical","s","Color",0.5*col,"MarkerFaceColor",0.5*col ,"MarkerSize",mark_size,"LineWidth",line_size)
    end

    errorbar(2 , mean_, std_+err_,"vertical","s","Color",0.8*col,"MarkerFaceColor",0.8*col ,"MarkerSize",mark_size,"LineWidth",line_size)
    errorbar(2 , mean_, std_,"vertical","s","Color",col,"MarkerFaceColor",col ,"MarkerSize",mark_size,"LineWidth",line_size)
    
    % soaked
    mean_ = mean(YoungsModuli(j,3:4));
    std_ = std(YoungsModuli(j,3:4));
    err_ = (1/2)*sqrt(sum(err_E(j,3:4).^2));

    col = 0.5 + 0.5*rand(3,1);
    for i=1:2 %size(samp_inds,2)
        errorbar(3 - w_space + 2*w_space*((i-1)/(2-1)), YoungsModuli(j,2+i), err_E(j,2+i),"vertical","s","Color",0.5*col,"MarkerFaceColor",0.5*col ,"MarkerSize",mark_size,"LineWidth",line_size)
    end

    errorbar(3 , mean_, std_+err_,"vertical","s","Color",0.8*col,"MarkerFaceColor",0.8*col ,"MarkerSize",mark_size,"LineWidth",line_size)
    errorbar(3 , mean_, std_,"vertical","s","Color",col,"MarkerFaceColor",col ,"MarkerSize",mark_size,"LineWidth",line_size)
    

    errorbar(0,Youngs_paper(j),err_Youngs_paper(j),"vertical","sk","MarkerFaceColor","k","MarkerSize",mark_size,"LineWidth",line_size)
    errorbar(1,Youngs_company(j),err_Youngs_company(j),"vertical","sk","MarkerFaceColor","k","MarkerSize",mark_size,"LineWidth",line_size)
    

    xlim([-1,4])
    xticks(0:3)
    xlabel(resins(j))
    ylabel("Young's Modulus [MPa]")
    xticklabels(["MTS value","Company Value","Control Beam","Soaked Beam"])
    xtickangle(22.5)
    ylim([0,3000])

end

if to_save
    saveas(gcf,"SummaryComparison.svg","svg")
end

function err = ModelError(E,x0,y0,theta0,params,S_data,m1,m2,x1,y1,x2,y2)

g = 9.81;
[x1_,y1_] = HingeDisplacement(E,x0,y0,theta0,params,-m1*g,S_data);
[x2_,y2_] = HingeDisplacement(E,x0,y0,theta0,params,-m2*g,S_data);

dx_err = 1000*[x1_ - x1; x2_ - x2];
dy_err = 1000*[y1_ - y1; y2_ - y2];
err = sum( dx_err.^2 + dy_err.^2 );

end