function littlerock


   %plot a sounding 
    filename='littlerock.nc';
    fprintf('reading file: %s\n',filename);
    file_struct=nc_info(filename);
    c=constants;
    %
    % grab the March 2 12Z sounding
    %
    sound_var = file_struct.Dataset(4).Name;
    fprintf('found sounding: %s\n',sound_var);
    press=nc_varget(filename,sound_var,[0,0],[Inf,1]);
    temp=nc_varget(filename,sound_var,[0,2],[Inf,1]);
    height = nc_varget(filename,sound_var,[0,1],[Inf,1]);
    dewpoint=nc_varget(filename,sound_var,[0,3],[Inf,1]);
    fh=figure(1);
    semilogy(temp,press);
    hold on;
    semilogy(dewpoint,press);
    set(gca,'yscale','log','ydir','reverse');
    ylim([400,1000]);
    ylabel('press (hPa)');
    xlabel('Temp (deg C)');
    title('sounding 1');
    hold off;
    figHandle=figure(2);
    skew=30.;
    [figHandle,outputws,handlews]=makeSkew(figHandle,skew);
    xtemp=convertTempToSkew(temp,press,skew);
    xdew=convertTempToSkew(dewpoint,press,skew);
    semilogy(xtemp,press,'g-','linewidth',5);
    semilogy(xdew,press,'b-','linewidth',5);
    %use lowest sounding level for adiabat
    thetaeVal=thetaes(temp(1) + c.Tc,press(1)*100.);
    [pressVals,tempVals]=calcAdiabat(press(1)*100.,thetaeVal,400.e2);
    xTemp=convertTempToSkew(tempVals - c.Tc,pressVals*1.e-2,skew);
    semilogy(xTemp,pressVals*1.e-2,'r-','linewidth',5);
    ylim([400,1000.]);
    xleft=convertTempToSkew(-20,1.e3,skew);
    xright=convertTempToSkew(25.,1.e3,skew);
    xlim([xleft,xright]);
    %
    % interpolator fails if two pressure values
    % are the same -- nudge them
    %
    newPress=nudgepress(press);
    interpTenv=@(pVals) interp1(newPress,temp,pVals);
    interpTdEnv=@(pVals) interp1(newPress,dewpoint,pVals);
    
    interpTenv_hgt=@(hVals) interp1(height,temp,hVals);
    interpPress_hgt=@(hVals) interp1(height,press,hVals);
    interpTdEnv_hgt=@(hVals) interp1(height,dewpoint,hVals);
    
    
    trytemp=interpTenv(pressVals*1.e-2);
    xTemp=convertTempToSkew(trytemp,pressVals*1.e-2,skew);
    semilogy(xTemp,pressVals*1.e-2,'b.','markersize',5);
    hold off;
    pressLevs=linspace(400,press(1),100)*1.e2;
    %
    % starting integrating from first sounding level
    %
    pressLevs=fliplr(pressLevs);
    TvDiffHandle=@(pVals) calcTvDiff(pVals,thetaeVal,interpTenv,interpTdEnv);
    for i=1:numel(pressLevs)
       Tvdiff(i)=TvDiffHandle(pressLevs(i));
    end
    
    
    
    
    %TvDiffHandle_B=@(hVals) B(hVals,press,thetae0,interpTenv,interpTdEnv,interpTenv,interpTdEnv);
      
    %for i=1:numel(pressLevs)
    %   B(i)=TvDiffHandle_B(pressLevs(i));
    %end
    
    
    figure(3);
    clf;
    plot(Tvdiff,pressLevs*0.01,'k-');
    set(gca,'ydir','reverse');
    ylabel('pressure (hPa)');
    xlabel('Virtual temperature difference (K)');
    title('Tvdiff vs. pressure');
    cumCAPE= -c.Rd*cumsum(Tvdiff(2:end).*diff(log(pressLevs)));
    figure(4);
    clf;
    plot(cumCAPE,pressLevs(2:end)*0.01,'k-');
    title('cumulative CAPE (J/kg) vs. pressure (hPa)');
    set(gca,'ydir','reverse');
    figure(5);
    clf;
    %
    % equate kinetic and potential energy toget maximum
    % updraft speed
    %
    maxvel=sqrt(2.*cumCAPE);
    plot(maxvel,pressLevs(2:end)*0.01,'k-');
    title('maximum updraft (m/s) vs. pressure (hPa)');
    set(gca,'ydir','reverse');
    
    %derives = @(t) 
    tspan=0:0.1:2500;
    yinit=[0,0];
    
    
   
    F = @(t,y) derives(t,y,press(1),thetaeVal,interpTenv_hgt,interpTdEnv_hgt);
    
    %Not right yet :(
    [t,y] = ode45(F, tspan, yinit);
   
   
    
    
    
end    
   
function newPress=nudgepress(pressVec)
    %if two balloon pressure levels are idential
    %add a factor of 0.1% to the second one
    %so interpolation will work
    newPress=pressVec;
    hit=find(abs(diff(newPress)) < 1.e-8);
    newPress(hit+1)=pressVec(hit) + 1.e-3*pressVec(hit);
end

function TvDiff=calcTvDiff(press,thetae0,interpTenv,interpTdEnv)
    %calcTvDiff(press,thetae0,interpTenv,interpTdenv)
    %input: press (Pa), thetae0 (K), plus function handles for T,Td soundings
    %output: TvDiff (K)
    %neglect liquid water loading in the virtual temperature
    c=constants;
    Tcloud=findTmoist(thetae0,press);
    wvcloud=wsat(Tcloud,press);
    Tvcloud=Tcloud*(1. + c.eps*wvcloud);
    Tenv=interpTenv(press*1.e-2) + c.Tc;
    Tdenv=interpTdEnv(press*1.e-2) + c.Tc;
    wvenv=wsat(Tdenv,press);
    Tvenv=Tenv*(1. + c.eps*wvenv);
    TvDiff=Tvcloud - Tvenv;
end

function buoy_out=B(press,thetaeVal,interpTenv,interpTdEnv)
    %find Buoyance, press here is the interp press
    c=constants;
    Tcloud=findTmoist(thetaeVal,press);
    wvcloud=wsat(Tcloud,press);
    Tvcloud=Tcloud*(1. + c.eps*wvcloud);
    Tenv=interpTenv(press*1.e-2) + c.Tc;
    Tdenv=interpTdEnv(press*1.e-2) + c.Tc;
    wvenv=wsat(Tdenv,press);
    Tvenv=Tenv*(1. + c.eps*wvenv);
    g=9.8;
    buoy_out=g*((Tvcloud - Tvenv)/Tvenv);
end

function yp=derives(t,y,press,thetaeVal,interpTenv,interpTdEnv)
  yp=zeros(2,1); % since output must be a column vector
  
  yp(1)=y(2);
  
     
  yp(2)= B(press,thetaeVal,interpTenv,interpTdEnv);
end

