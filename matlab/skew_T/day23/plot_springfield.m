function plot_springfield()

%plot a sounding from soundings.nc
clear all
close all
 filename='springfield.nc';
 file_struct=nc_info(filename);
 c=constants;
 %
 % grap the first sounding pressure and temperature
 %
 sound_var = file_struct.Dataset(3).Name;
 press=nc_varget(filename,sound_var,[0,0],[Inf,1]);
 temp=nc_varget(filename,sound_var,[0,2],[Inf,1]);
 dewpoint=nc_varget(filename,sound_var,[0,3],[Inf,1]);
 fh=figure(1);
 semilogy(temp,press)
 hold on;
 semilogy(dewpoint,press)
 set(gca,'yscale','log','ydir','reverse');
 ylim([400,1000]);
 ylabel('press (hPa)')
 xlabel('Temp (deg C)')
 title('sounding 1')
 hold off;
 figHandle=figure(2);
  skew=30.;
 [figHandle,outputws,handlews]=makeSkew(figHandle,skew);
 xtemp=convertTempToSkew(temp,press,skew);    
 xdew=convertTempToSkew(dewpoint,press,skew);    
 semilogy(xtemp,press,'g-','linewidth',5);
 semilogy(xdew,press,'b-','linewidth',5);
 [xTemp,thePress]=ginput(1);
 Tclick=convertSkewToTemp(xTemp,thePress,skew);    
 thetaeVal=thetaes(Tclick + c.Tc,thePress*100.);
 fprintf('ready to draw moist adiabat, thetae=%8.2f\n',thetaeVal);
 ylim([400,1000.])
 
 %plotting theta that was picked by ginput
 
 w=wsat(Tclick+c.Tc, thePress*100);
 
 press_plot=thePress:-10:300;
 for i = 1:length(press_plot)
     
    [temp_thetae(i), wl(i), wv(i)]=tinvert_thetae(thetaeVal,w,press_plot(i)*100);
 end
 thetaeplot2 = convertTempToSkew(temp_thetae-c.Tc, press_plot, skew);
 a1=plot(thetaeplot2, press_plot, 'k-','linewidth',3);
 
 
% plotting theta for LCL and CAPE calc
 [Tlcl_surf, plcl_surf] = LCLfind(dewpoint(1)+c.Tc, temp(1)+c.Tc, press(1)*100);
 w=wsat(Tlcl_surf+c.Tc, plcl_surf*100);
 thetae_lcl = thetaes(Tlcl_surf, plcl_surf);
 press_plot=1000:-10:300;
 for i = 1:length(press_plot)
     
    [temp_thetae(i), wl(i), wv(i)]=tinvert_thetae(thetae_lcl,w,press_plot(i)*100);
 end
 thetaeplot2 = convertTempToSkew(temp_thetae-c.Tc, press_plot, skew);
 a2=plot(thetaeplot2, press_plot, 'r-','linewidth',5);
 legend(a1,'Thetae clicked')
 legend([a1 a2],'Thetae clicked','Thetae LCL')
 
 
 thetae_750 = temp_thetae(find(press_plot==750));%in K
 temp_750 = temp(find(press>740&press<759))+c.Tc;%in K
  
tmp_750=[thetae_750:0.0325:temp_750];%not skewed yet
 for i=1:length(tmp_750)
     press_btm(i) = 750;
 end
 
 thetae_400 = temp_thetae(find(press_plot==400));%in K
 temp_400 = temp(find(press==400))+c.Tc;%in K
tmp_400=[temp_400:.125:thetae_400];%not skewed yet
 for i=1:length(tmp_400)
     press_top(i) = 400;
 end
 temp(27)=temp(27)+0.1;
 
 theta_area = [tmp_400(1:end-1) temp_thetae(find(press_plot==400):-1:find(press_plot==750)) tmp_750(2:end)];
 temp_area = temp(find(press==400):-1:find(press>740&press<759))+c.Tc;
 temp_area = temp_area';
 theta_area(end)=temp_area(end);
 press(find(press>740&press<759))=750;
 temp(find(press>740&press<759))=tinvert_thetae(temp_thetae(find(press_plot==750)), w, 750*100);
 press_temp = press(find(press==400):-1:find(press>740&press<759))';
 press_theta = [press_top(1:end-1) press_plot(find(press_plot==400):-1:(find(press_plot==750))) press_btm(2:end)];
 
 CAPE=plot_work(theta_area, press_theta*100, temp_area, press_temp*100)
 
 %press_theta = 
 %plot(tmp, press_btm)
 %calculate area:
 %cape=plot_work(
 
 
 hold off;
 

% $$$         double Mar-17-2011-00Z(dim_138, var_cols) ;
% $$$         double Mar-17-2011-12Z(dim_139, var_cols) ;
% $$$         double Mar-18-2011-00Z(dim_128, var_cols) ;
% $$$         double Mar-18-2011-12Z(dim_142, var_cols) ;
% $$$         double Mar-19-2011-00Z(dim_39, var_cols) ;

end
function netWork=plot_work(alphaCAB,pressCAB,alphaBC,pressBC)
    figure(3);
    %clf;
    interpCABfun=@(alphavals) interp1(alphaCAB,pressCAB,alphavals);
    interpBCfun=@(alphavals) interp1(alphaBC,pressBC,alphavals);
    interpAlphaVals=linspace(247,286);
    interpCAB=interpCABfun(interpAlphaVals);
    interpBC=interpBCfun(interpAlphaVals);
    plot(alphaCAB,pressCAB*1.e-2);
    ylim([400,1000.]);
    hold on;
    plot(alphaBC,pressBC*1.e-2);
    plot(interpAlphaVals,interpCAB*1.e-2,'r*');
    plot(interpAlphaVals,interpBC*1.e-2,'g*');
    start=alphaBC(end);
    stop=alphaBC(1);
    workBC=quad(interpBCfun,start,stop);
    workCAB=quad(interpCABfun,start, stop);
    netWork=abs(workCAB) - abs(workBC);
    out_mesg={'workBC %9.3g (J/kg)\n',...
          'workCAB %9.3g (J/kg)\n',...
          'net work done by system %9.3g (J/kg)= CAPE\n'};
    fprintf(strcat(out_mesg{:}),workBC,workCAB,abs(workCAB) - abs(workBC));
    title('pressure - temperature plot');
    hold off;
end        