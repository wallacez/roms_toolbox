function varargout = compare_gui(varargin)
% COMPARE_GUI MATLAB code for compare_gui.fig
%      COMPARE_GUI, by itself, creates a new COMPARE_GUI or raises the existing
%      singleton*.
%
%      H = COMPARE_GUI returns the handle to a new COMPARE_GUI or the handle to
%      the existing singleton*.
%
%      COMPARE_GUI('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in COMPARE_GUI.M with the given input arguments.
%
%      COMPARE_GUI('Property','Value',...) creates a new COMPARE_GUI or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before compare_gui_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to compare_gui_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help compare_gui

% Last Modified by GUIDE v2.5 19-Mar-2019 17:13:23

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @compare_gui_OpeningFcn, ...
                   'gui_OutputFcn',  @compare_gui_OutputFcn, ...
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


% --- Executes just before compare_gui is made visible.
function compare_gui_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to compare_gui (see VARARGIN)

% Choose default command line output for compare_gui
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% Clear all variables (including globals)
clear all

% UIWAIT makes compare_gui wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = compare_gui_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


% --- Executes on selection change in listbox_exp.
function listbox_exp_Callback(hObject, eventdata, handles)
% hObject    handle to listbox_exp (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns listbox_exp contents as cell array
%        contents{get(hObject,'Value')} returns selected item from listbox_exp
    exp_selected = hObject.String(hObject.Value);
    
    if strcmp(exp_selected,'Reset')
    % Repopulate listbox with all available
    % experiments.
        listbox_exp_CreateFcn(hObject,eventdata,handles)
    else
    % Resize the size of hObject so it only holds reset and the selected
    % experiment.
        val=[hObject.String(1);exp_selected];
        hObject.Value=2;
        hObject.String = val;
      
        % Get list of files in selected experimental directory
        flist = get_flist_averages(exp_selected);
        
        % Get lat, lon, and vertical levels info
        M = get_grid_vars();  %<--- Pass this container around to various methods
        
        % Get dates from selected experiment and populate the "Date"
        % listbox.
        ot = combine_var_times(flist, 'ocean_time',0);
        if length(ot) == 1
            % Only 1 timefile read in
            ot = datestr(ot./86400,1);
        else
            ot = datestr(ot(2:end)./86400,1);  % Format as dd-mm-yyyy 
        end
        handles.listbox_date.String = ot;
        
        % set slider endpoints
        handles.slider_play.Min = 1;
        handles.slider_play.Max = length(ot);
        handles.slider_play.Value = handles.slider_play.Min;
        
        % Generate list of variables available for plotting.
        % Get only relevant physical and biological variables, which
        % usually is everything after the 'ocean_time' variable.
        nc = ncinfo(flist{1});
        varnamelist = {nc.Variables.Name};
        
        % Add combined variables, e.g. log10 Chls + ChlD
        idx_end = length(varnamelist);
        varnamelist{1,idx_end+1} = 'Log10 Total Chl';
        varnamelist{1,idx_end+2} = 'Log10 Total Chl 1.59';
        varnamelist{1,idx_end+3} = 'Total Chl';
        varnamelist{1,idx_end+4} = 'Log10 Total Phy';
        varnamelist{1,idx_end+5} = 'Total Phy';
        varnamelist{1,idx_end+6} = 'Log10 ChlS';
        varnamelist{1,idx_end+7} = 'Log10 ChlS 1.59';
        varnamelist{1,idx_end+8} = 'Log10 ChlD';
        varnamelist{1,idx_end+9} = 'Log10 ChlD 1.59';
        varnamelist{1,idx_end+10} = 'Chl:N Tot';
        varnamelist{1,idx_end+11} = 'Chl:N S';
        varnamelist{1,idx_end+12} = 'Chl:N L';    
        
        idx = find(strcmp(varnamelist,'ocean_time'));
        
        % Populate popup menus with variable lists.
        % Menu 2 will be invisible unless "2 plots" button is selected.
        handles.popupmenu_var1.String = varnamelist(idx+1:end);
        handles.popupmenu_var2.String = varnamelist(idx+1:end);
        
     end
    

% --- Executes during object creation, after setting all properties.
function listbox_exp_CreateFcn(hObject, eventdata, handles)
% hObject    handle to listbox_exp (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: listbox controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
% Populate experiment list with directories on morrison
    d=dir('/home/morrison/data2/zach');
    dirs=find([d.isdir]==1);  % find all in directory that are only folders
    exp_list={d(dirs(2:end)).name};  % exclude '.' directory
    exp_list{1}='Reset';  % Include option to choose a different experiment
    hObject.String = exp_list;


% --- Executes on selection change in listbox_date.
function listbox_date_Callback(hObject, eventdata, handles)
% hObject    handle to listbox_date (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns listbox_date contents as cell array
%        contents{get(hObject,'Value')} returns selected item from
%        listbox_date

    
% --- Executes during object creation, after setting all properties.
function listbox_date_CreateFcn(hObject, eventdata, handles)
% hObject    handle to listbox_date (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: listbox controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end    


        % --- Executes on selection change in popupmenu_var1.
function popupmenu_var1_Callback(hObject, eventdata, handles)
% hObject    handle to popupmenu_var1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns popupmenu_var1 contents as cell array
%        contents{get(hObject,'Value')} returns selected item from popupmenu_var1
    global VARDATA1 VAR1_SELECTED h VAR1_DEPTH Z_REAL X_LOC Y_LOC
    exp_dir = handles.listbox_exp.String(2);
    flist = get_flist_averages(exp_dir);
    
    % Read selected variable and build time series at surface
    VAR1_SELECTED = hObject.String{hObject.Value};
    VARDATA1 = combine_var_times(flist, VAR1_SELECTED,0);
    
    % Obtain selected date
    h_dates = handles.listbox_date;
    curr_date = h_dates.Value;  % index of user-selected date
    num_dates = length(h_dates.String(:,1));  % all dates in experiment
    
    % Set slider endpoints.
    h_slider = handles.slider_play;
    h_slider.Min = 1;
    h_slider.Max = num_dates;
    
    if strncmp(VAR1_SELECTED,'Log10',5)
        VARDATA1(VARDATA1<0)=NaN;  % temporary fix
    end
    
    % Plot selected date
    h_plot2button = handles.radiobutton_plot2;
    if h_plot2button(1).Value == 0 && h_plot2button(2).Value == 1
        % 2 horizontal plots
        ax = handles.axes_2horz1;
        cla(ax,'reset');colorbar(ax,'delete');  % reset axes
        
        VAR1_SELECTED
        if strncmp(VAR1_SELECTED,'Log10',5)
            % Log10 Transform 
            pcolor(ax,log10(VARDATA1(:,:,curr_date)'));
            shading(ax,'flat');colorbar('peer',ax);  % Need to specify axes here.
            caxis(ax,[-1 0.7]);
        else
            pcolor(ax,VARDATA1(:,:,curr_date)');
            shading(ax,'flat');colorbar('peer',ax);  % Need to specify axes here
            caxis(ax,'auto');
        end
    else
        % 1 horizontal plot
        ax = handles.axes_1horz;
        cla(ax,'reset');colorbar(ax,'delete');  % reset axes
        
        VAR1_SELECTED
        if strncmp(VAR1_SELECTED,'Log10',5)
            % Log10 Transform 
            pcolor(ax,log10(VARDATA1(:,:,curr_date)'));
            shading(ax,'flat');colorbar('peer',ax);
            caxis(ax,[-1 0.7]);
        else
            pcolor(ax,VARDATA1(:,:,curr_date)');
            shading(ax,'flat');colorbar('peer',ax);  
        end
    end

    % Set color scale
    handles.edit_caxis_plot1_min.String = round(ax.CLim(1));        
    if strncmp(VAR1_SELECTED,'Log10',5)
        handles.edit_caxis_plot1_max.String = round(ax.CLim(2),1);
    else
        handles.edit_caxis_plot1_max.String = round(ax.CLim(2));
    end    
    colormap(ax,jet);
    
    % Plot 200-m isobath.  Fix colorbar limits so they are not overwritten
    % by the contour plot.
    caxis(ax,[ax.CLim(1) ax.CLim(2)]);
    hold(ax,'on');contour(ax,h',[200 200],'w','Linewidth',2);
    title(ax,h_dates.String(curr_date,:))

    % Change vertical plot associated with variable 1, if section has been
    % chosen.
    v_plotbutton = handles.radiobutton_vertPlotY.Value;    
    if v_plotbutton == 1 && ~isempty(VAR1_DEPTH)    
        % Read in new depth data
        VAR1_DEPTH = combine_var_times(flist,VAR1_SELECTED,1);
        if strncmp(VAR1_SELECTED,'Log10',5)
            VAR1_DEPTH(VAR1_DEPTH<0)=NaN;
        end        
        
        if h_plot2button(1).Value == 0 && h_plot2button(2).Value == 1
            % Plot to left axes since 2 depth variables to plot
            ax1 = handles.axes_2vert1;
        else
            % Only single depth plot
            ax1 = handles.axes_1vert;
        end
        if strncmp(VAR1_SELECTED,'Log10',5)
            % Log10 Transform 
            pcolor(ax1, ...
                       repmat(X_LOC:482,40,1), ...
                       squeeze(Z_REAL(X_LOC:482,Y_LOC,:))', ...
                       squeeze(log10(VAR1_DEPTH(:,:,curr_date)')));
            shading(ax1,'interp');colorbar('peer',ax1,'Location','Southoutside');  % Need to specify axes here.
            caxis(ax1,[-1 0.7]);
        else
            pcolor(ax1, ...
                       repmat(X_LOC:482,40,1), ...
                       squeeze(Z_REAL(X_LOC:482,Y_LOC,:))', ...
                       squeeze(VAR1_DEPTH(:,:,curr_date)'));
            shading(ax1,'interp');colorbar('peer',ax1,'Location','Southoutside');  % Need to specify axes here
            caxis(ax1,'auto');           
        end
        
        % default to currently specified depth
        ax1.YLim = [str2double(handles.edit_depthMax.String) ...
                            str2double(handles.edit_depthMin.String)];
        
        % Set color scale
        handles.edit_depth1ColorbarMin.String = round(ax1.CLim(1));        
        if strncmp(VAR1_SELECTED,'Log10',5)
            handles.edit_depth1ColorbarMax.String = round(ax1.CLim(2),1);
        else
            handles.edit_depth1ColorbarMax.String = round(ax1.CLim(2));
        end
        colormap(ax1,jet);        
    end    


% --- Executes during object creation, after setting all properties.
function popupmenu_var1_CreateFcn(hObject, eventdata, handles)
% hObject    handle to popupmenu_var1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in popupmenu_var2.
function popupmenu_var2_Callback(hObject, eventdata, handles)
% hObject    handle to popupmenu_var2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns popupmenu_var2 contents as cell array
%        contents{get(hObject,'Value')} returns selected item from popupmenu_var2
    global VARDATA2 VAR2_SELECTED VAR2_DEPTH Z_REAL X_LOC Y_LOC h
    exp_dir = handles.listbox_exp.String(2);
    flist = get_flist_averages(exp_dir);
    
    % Read selected variable and build time series at surface
    VAR2_SELECTED = hObject.String{hObject.Value};
    VARDATA2 = combine_var_times(flist, VAR2_SELECTED,0);
    
    % Obtain selected date
    h_dates = handles.listbox_date;
    curr_date = h_dates.Value;  % index of user-selected date
    
    % Plot selected date
    ax = handles.axes_2horz2;
    cla(ax,'reset');colorbar(ax,'delete');  % reset axes
    
    VAR2_SELECTED
    
     if strncmp(VAR2_SELECTED,'Log10',5)
        VARDATA2(VARDATA2<0)=NaN;  % temporary fix
    end
    
    if strncmp(VAR2_SELECTED,'Log10',5)
        % Log10 Transform 
        pcolor(ax,log10(VARDATA2(:,:,curr_date)'));
        shading(ax,'flat');colorbar('peer',ax);  % Need to specify axes, here idk why.
        caxis(ax,[-1 0.7]);
    else
        pcolor(ax,VARDATA2(:,:,curr_date)');
        shading(ax,'flat');colorbar('peer',ax);  % Need to specify axes, here idk why.
    end

    % Set color scale
    handles.edit_caxis_plot2_min.String = round(ax.CLim(1));        
    if strncmp(VAR2_SELECTED,'Log10',5)
        handles.edit_caxis_plot2_max.String = round(ax.CLim(2),1);
    else
        handles.edit_caxis_plot2_max.String = round(ax.CLim(2));
    end     
    colormap(ax,jet);
    
    % Plot 200-m isobath.  Set colorbar limits so they are not overwritten
    % by the contour plot.
    caxis(ax,[ax.CLim(1) ax.CLim(2)]);
    hold(ax,'on');contour(ax,h',[200 200],'w','Linewidth',2);
    title(ax,h_dates.String(curr_date,:))    
    
    % Change vertical plot associated with variable 2, if section has been
    % chosen.
    v_plotbutton = handles.radiobutton_vertPlotY.Value;
    if v_plotbutton == 1 && ~isempty(VAR2_DEPTH)        
        % Read in new depth data
        VAR2_DEPTH = combine_var_times(flist,VAR2_SELECTED,1);
        if strncmp(VAR2_SELECTED,'Log10',5)
            VAR2_DEPTH(VAR2_DEPTH<0)=NaN;
        end      
        
        % Plot second variable
        ax2 = handles.axes_2vert2;
        if strncmp(VAR2_SELECTED,'Log10',5)
            % Log10 Transform 
            pcolor(ax2, ...
                       repmat(X_LOC:482,40,1), ...
                       squeeze(Z_REAL(X_LOC:482,Y_LOC,:))', ...
                       squeeze(log10(VAR2_DEPTH(:,:,curr_date)')));
            shading(ax2,'interp');colorbar('peer',ax2,'Location','Southoutside');  % Need to specify axes here.
            caxis(ax2,[-1 0.7]);
        else
            pcolor(ax2, ...
                       repmat(X_LOC:482,40,1), ...
                       squeeze(Z_REAL(X_LOC:482,Y_LOC,:))', ...
                       squeeze(VAR2_DEPTH(:,:,curr_date)'));
            shading(ax2,'interp');colorbar('peer',ax2,'Location','Southoutside');  % Need to specify axes here
            caxis(ax2,'auto')
        end
        
        % default to currently specified depth
        ax2.YLim = [str2double(handles.edit_depthMax.String) ...
                            str2double(handles.edit_depthMin.String)];
        ax2.YTickLabels = '';  % Turn off for clarity
        
        % Set color scale
        handles.edit_depth2ColorbarMin.String = round(ax2.CLim(1));        
        if strncmp(VAR2_SELECTED,'Log10',5)
            handles.edit_depth2ColorbarMax.String = round(ax2.CLim(2),1);
        else
            handles.edit_depth2ColorbarMax.String = round(ax2.CLim(2));
        end        
        colormap(ax2,jet);        
    end


% --- Executes during object creation, after setting all properties.
function popupmenu_var2_CreateFcn(hObject, eventdata, handles)
% hObject    handle to popupmenu_var2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
    hObject.Visible = 'Off';  % default to only plotting 1 variable.
    


% --- Executes on slider movement.
function slider_play_Callback(hObject, eventdata, handles)
% hObject    handle to slider_play (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider
   
    global VARDATA1 VAR1_SELECTED VARDATA2 VAR2_SELECTED ...
               VAR1_DEPTH VAR2_DEPTH Z_REAL X_LOC Y_LOC h
    
    h_plot2button = handles.radiobutton_plot2;
    v_plotbutton = handles.radiobutton_vertPlotY.Value;
    
    % set slider endpoints
    hObject.Min = 1;
    hObject.Max = length(handles.listbox_date.String(:,1));
    curr_date = round(hObject.Value);
    
    % set date in date_box to current slider value
    handles.listbox_date.Value = curr_date;
    
    % get caxis handles
    cmin_h1 = str2double(handles.edit_caxis_plot1_min.String);
    cmax_h1 = str2double(handles.edit_caxis_plot1_max.String);
    if v_plotbutton == 1
        cmin_v1 = str2double(handles.edit_depth1ColorbarMin.String);
        cmax_v1 = str2double(handles.edit_depth1ColorbarMax.String);
    end
    if h_plot2button(1).Value == 0 && h_plot2button(2).Value == 1
        cmin_h2 = str2double(handles.edit_caxis_plot2_min.String);
        cmax_h2 = str2double(handles.edit_caxis_plot2_max.String);
        if v_plotbutton == 1
            cmin_v2 = str2double(handles.edit_depth2ColorbarMin.String);
            cmax_v2 = str2double(handles.edit_depth2ColorbarMax.String);
        end
    end
    
    % Plot data at appropriate slider point
    if h_plot2button(1).Value == 0 && h_plot2button(2).Value == 1
        % 2 horizontal plots
        ax_h1 = handles.axes_2horz1;
        ax_h2 = handles.axes_2horz2;
        
        % Prevent slider value from exceeding the date range of the 
        % shortest run being plotted.
        size_var1 = size(VARDATA1);
        size_var2 = size(VARDATA2);
        if curr_date > size_var1(end)
            curr_date = size_var1(end);
        elseif curr_date > size_var2(end)
            curr_date = size_var2(end);
        end
        
        % Plot first variable
        cla(ax_h1,'reset');colorbar(ax_h1,'delete');  % reset axes
        if strncmp(VAR1_SELECTED,'Log10',5)
            % Log10 Transform 
            pcolor(ax_h1,log10(VARDATA1(:,:,curr_date)'));
        else 
            pcolor(ax_h1,VARDATA1(:,:,curr_date)');
        end
        shading(ax_h1,'flat');colorbar('peer',ax_h1);
        caxis(ax_h1,[cmin_h1 cmax_h1]);
        hold(ax_h1,'on');contour(ax_h1,h',[200 200],'w','Linewidth',2);
        title(ax_h1,handles.listbox_date.String(curr_date,:))        
        colormap(ax_h1,jet);
        
        % Plot second variable
        cla(ax_h2,'reset');colorbar(ax_h2,'delete');  % reset axes
        if strncmp(VAR2_SELECTED,'Log10',5)
            % Log10 Transform 
            pcolor(ax_h2,log10(VARDATA2(:,:,curr_date)'));
        else
            pcolor(ax_h2,VARDATA2(:,:,curr_date)');
        end
        shading(ax_h2,'flat');colorbar('peer',ax_h2);
        caxis(ax_h2,[cmin_h2 cmax_h2]);
        hold(ax_h2,'on');contour(ax_h2,h',[200 200],'w','Linewidth',2);
        title(ax_h2,handles.listbox_date.String(curr_date,:));        
        colormap(ax_h2,jet);
       
        % Plot first vertical variable, if appropriate
        if v_plotbutton == 1 && ~isempty(VAR1_DEPTH) && ~isempty(VAR2_DEPTH)
            % Plot to left axes since 2 depth variables to plot
            ax1 = handles.axes_2vert1;
            if strncmp(VAR1_SELECTED,'Log10',5)
                % Log10 Transform 
                pcolor(ax1, ...
                           repmat(X_LOC:482,40,1), ...
                           squeeze(Z_REAL(X_LOC:482,Y_LOC,:))', ...
                           squeeze(log10(VAR1_DEPTH(:,:,curr_date)')));
                shading(ax1,'interp');colorbar('peer',ax1,'Location','Southoutside');  % Need to specify axes here.
                caxis(ax1,[-1 0.7]);
            else
                pcolor(ax1, ...
                           repmat(X_LOC:482,40,1), ...
                           squeeze(Z_REAL(X_LOC:482,Y_LOC,:))', ...
                           squeeze(VAR1_DEPTH(:,:,curr_date)'));
                shading(ax1,'interp');colorbar('peer',ax1,'Location','Southoutside');  % Need to specify axes here
                caxis(ax1,[cmin_v1 cmax_v1]);
            end
            % default to currently specified depth
            ax1.YLim = [str2double(handles.edit_depthMax.String) ...
                                str2double(handles.edit_depthMin.String)];
            colormap(ax1,jet);                 
        
        % Plot second vertical variable       
            ax2 = handles.axes_2vert2;
            if strncmp(VAR2_SELECTED,'Log10',5)
                % Log10 Transform 
                pcolor(ax2, ...
                           repmat(X_LOC:482,40,1), ...
                           squeeze(Z_REAL(X_LOC:482,Y_LOC,:))', ...
                           squeeze(log10(VAR2_DEPTH(:,:,curr_date)')));
                shading(ax2,'interp');colorbar('peer',ax2,'Location','Southoutside');  % Need to specify axes here.
                caxis(ax2,[-1 0.7]);
            else
                pcolor(ax2, ...
                           repmat(X_LOC:482,40,1), ...
                           squeeze(Z_REAL(X_LOC:482,Y_LOC,:))', ...
                           squeeze(VAR2_DEPTH(:,:,curr_date)'));
                shading(ax2,'interp');colorbar('peer',ax2,'Location','Southoutside');  % Need to specify axes here
                caxis(ax2,[cmin_v2 cmax_v2]);
            end
            % default to currently specified depth
            ax2.YLim = [str2double(handles.edit_depthMax.String) ...
                                str2double(handles.edit_depthMin.String)];
            ax2.YTickLabels = '';  % Turn off for clarity
            colormap(ax2,jet); 
        end
    else
        % Only 1 plot
        ax = handles.axes_1horz;
        if strncmp(VAR1_SELECTED,'Log10',5)
            % Log10 Transform 
            pcolor(ax,log10(VARDATA1(:,:,curr_date)'));shading flat;colorbar;
            title(handles.listbox_date.String(curr_date,:));
        else
            pcolor(ax,VARDATA1(:,:,curr_date)');shading flat;colorbar;
            caxis([cmin_h1 cmax_h1]);
            title(handles.listbox_date.String(curr_date,:))
        end
        colormap(ax,jet);        
    end   
    

% --- Executes during object creation, after setting all properties.
function slider_play_CreateFcn(hObject, eventdata, handles)
% hObject    handle to slider_play (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end


% --- Executes on button press in togglebutton_play.
function togglebutton_play_Callback(hObject, eventdata, handles)
% hObject    handle to togglebutton_play (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

    h_plot2button = handles.radiobutton_plot2;
    
    % set variables containing handles to appropriate UI controls to be
    % passed to play_movie function.
    h_slider = handles.slider_play;
    h_dates = handles.listbox_date;
    h_cmin1 = handles.edit_caxis_plot1_min;
    h_cmax1 = handles.edit_caxis_plot1_max;
    
    if h_plot2button(1).Value == 0 && h_plot2button(2).Value == 1
        % 2 horizontal plots
        h_axes21 = handles.axes_2horz1;
        h_axes22 = handles.axes_2horz2;
        h_cmin2 = handles.edit_caxis_plot2_min;
        h_cmax2 = handles.edit_caxis_plot2_max;
    else % 1 horizontal plot
        h_axes1 = handles.axes_1horz;
    end

    if strcmp(hObject.String,'Resume')
        % Resume movie from where it was paused.
        curr_date = round(h_slider.Value);
    else
        % Start movie from user-selected date.
        curr_date = h_dates.Value;  % index of user-selected date
    end
    num_dates = length(h_dates.String(:,1));  % all dates in experiment
    
    % Set slider endpoints.
    h_slider.Min = 1;
    h_slider.Max = num_dates;
    
    % Toggle button string depending on depressed status.
    if hObject.Value == hObject.Max  % Button depressed, play movie activated
        hObject.String = 'Pause';
         if h_plot2button(1).Value == 0 && h_plot2button(2).Value == 1
             % 2 horizontal plots
             play_movie(h_slider, h_dates, curr_date, num_dates, ...
                                h_axes21, h_cmin1, h_cmax1,...
                                h_axes22, h_cmin2,h_cmax2)
         else
             % 1 horizontal plot
            play_movie(h_slider, h_dates, curr_date, num_dates, ...
                               h_axes1,h_cmin1,h_cmax1)
         end
        hObject.String = 'Play';
        hObject.Value = hObject.Min;
    else
        % Pause the animation
        hObject.String = 'Resume';
        pause
    end
    
    
function play_movie(h_slider, h_dates, start_date, end_date, varargin)
%
% Plot animation on provided axes (axes:  TODO)
% Input: h_slider(obj):  handle to slider_play UI control
%           h_dates(obj):  handle to listbox_dates UI control
%           start_date(int):  starting time index from listbox_date
%           end_date(int):    ending time indext from listbox_date
%           varargin:  Extra arguments account for variable number of axes
%                           passed to the function.
% Output: None
%     
    global VAR1_SELECTED VAR2_SELECTED VARDATA1 VARDATA2 h
    
    if nargin == 7  % 1 horizontal plot
        vartitle = h_dates.String;
        
        % Get colorscale range
        cmin = str2double(varargin{2}.String);
        cmax = str2double(varargin{3}.String);

        % Get axis to plot to
        ax_h1 = varargin{1};
        if strncmp(VAR1_SELECTED,'Log10',5)
            for i = start_date:end_date
                cla(ax_h1,'reset');colorbar(ax_h1,'delete');
                % update plot
                pcolor(ax_h1,log10(VARDATA1(1:480,1:768,i)'));shading flat;colorbar
                caxis([cmin cmax]); title(vartitle(i,:));
                hold(ax_h1,'on');contour(ax_h1,h',[200 200],'w','Linewidth',2);
                colormap(ax_h1,jet);
                % update slider position
                h_slider.Value = i;
                pause(1)
            end
        else
            for i = start_date:end_date
                cla(ax_h1,'reset');colorbar(ax_h1,'delete');
                % update plot
                pcolor(ax_h1,VARDATA1(1:480,1:768,i)');shading flat;colorbar;
                caxis([cmin cmax]); title(vartitle(i,:));
                hold(ax_h1,'on');contour(ax_h1,h',[200 200],'w','Linewidth',2);
                colormap(ax_h1,jet);
                % update slider position
                h_slider.Value = i;
                pause(1)
            end   
        end
        
    elseif nargin == 10 % 2 horizontal plots
        vartitle = h_dates.String;
        
        % Get colorscale range
        cmin1 = str2double(varargin{2}.String);
        cmax1 = str2double(varargin{3}.String);
        cmin2 = str2double(varargin{5}.String);
        cmax2 = str2double(varargin{6}.String);

        % Get axes to plot to
        ax_h1 = varargin{1};
        ax_h2 = varargin{4};
        
        % Prevent slider value from exceeding the
        %  date range of the shortest run being plotted.
        size_var1 = size(VARDATA1);
        size_var2 = size(VARDATA2);
        if end_date > size_var1(end)
            end_date = size_var1(end);
        elseif end_date > size_var2(end)
            end_date = size_var2(end);
        end
        
        % update plots
        for i = start_date:end_date
            % Plot 1
            cla(ax_h1,'reset');colorbar(ax_h1,'delete');
            if strncmp(VAR1_SELECTED,'Log10',5)
                pcolor(ax_h1,log10(VARDATA1(1:480,1:768,i)'));
            else
                pcolor(ax_h1,VARDATA1(1:480,1:768,i)');
            end
            shading(ax_h1,'flat');colorbar('peer',ax_h1);
            caxis(ax_h1,[cmin1 cmax1]); title(ax_h1,vartitle(i,:));
            hold(ax_h1,'on');contour(ax_h1,h',[200 200],'w','Linewidth',2);
            colormap(ax_h1,jet);
            
            % Plot 2
            cla(ax_h2,'reset');colorbar(ax_h2,'delete');
            if strncmp(VAR2_SELECTED,'Log10',5)
                pcolor(ax_h2,log10(VARDATA2(1:480,1:768,i)'));
            else
                pcolor(ax_h2,VARDATA2(1:480,1:768,i)');
            end
            shading(ax_h2,'flat');colorbar('peer',ax_h2);
            caxis(ax_h2,[cmin2 cmax2]); title(ax_h2,vartitle(i,:));
            hold(ax_h2,'on');contour(ax_h2,h',[200 200],'w','Linewidth',2);
            colormap(ax_h2,jet);
                
            % update slider position
            h_slider.Value = i;
            pause(1)
        end
    end
    
    
% --- Executes on button press in pushbutton_Section.
function pushbutton_Section_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_Section (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
    global X_LOC Y_LOC VAR1_SELECTED VAR2_SELECTED Z_REAL ...
               VAR1_DEPTH VAR2_DEPTH
    
    if ~exist('Z','var')
        z_r = load('/home/morrison/data2/zach/real_depths.mat','Z');
        Z_REAL = z_r.Z;
    end

    [X_LOC,Y_LOC] = ginput(1);
    X_LOC = round(X_LOC);
    Y_LOC = round(Y_LOC);
    handles.edit_gridPointX.String = int2str(X_LOC);
    handles.edit_gridPointY.String = int2str(Y_LOC);
    
    exp_dir = handles.listbox_exp.String(2);
    flist = get_flist_averages(exp_dir);
    
    % For checking if 1 or 2 vertical plots.
    h_plot2button = handles.radiobutton_plot2;
    
    % Load variable(s) with depths 
    VAR1_DEPTH = combine_var_times(flist, VAR1_SELECTED, 1);
    if h_plot2button(1).Value == 0 && h_plot2button(2).Value == 1
        % 2 vertical plots
        VAR2_DEPTH = combine_var_times(flist,VAR2_SELECTED,1);
    end
    
    % Check for negatives that would prevent pcolor plot of log10 transform
    if strncmp(VAR1_SELECTED,'Log10',5)
       VAR1_DEPTH(VAR1_DEPTH<0)=NaN;
    end
    if strncmp(VAR2_SELECTED,'Log10',5)
       VAR2_DEPTH(VAR2_DEPTH<0)=NaN;
    end        
    
    % Obtain selected date
    h_dates = handles.listbox_date;
    curr_date = h_dates.Value;  % index of user-selected date
    num_dates = length(h_dates.String(:,1));  % all dates in experiment
   
    % Plot selected date
    if h_plot2button(1).Value == 0 && h_plot2button(2).Value == 1
        % Plot first variable
        ax1 = handles.axes_2vert1;
        if strncmp(VAR1_SELECTED,'Log10',5)
            % Log10 Transform 
            pcolor(ax1, ...
                       repmat(X_LOC:482,40,1), ...
                       squeeze(Z_REAL(X_LOC:482,Y_LOC,:))', ...
                       squeeze(log10(VAR1_DEPTH(:,:,curr_date)')));
            shading(ax1,'interp');colorbar('peer',ax1,'Location','Southoutside');  % Need to specify axes here.
            caxis(ax1,[-1 0.7]);
        else
            pcolor(ax1, ...
                       repmat(X_LOC:482,40,1), ...
                       squeeze(Z_REAL(X_LOC:482,Y_LOC,:))', ...
                       squeeze(VAR1_DEPTH(:,:,curr_date)'));
            shading(ax1,'interp');colorbar('peer',ax1,'Location','Southoutside');  % Need to specify axes here
            caxis(ax1,'auto')
        end
        % Plot second variable
        ax2 = handles.axes_2vert2;
        if strncmp(VAR2_SELECTED,'Log10',5)
            % Log10 Transform 
            pcolor(ax2, ...
                       repmat(X_LOC:482,40,1), ...
                       squeeze(Z_REAL(X_LOC:482,Y_LOC,:))', ...
                       squeeze(log10(VAR2_DEPTH(:,:,curr_date)')));
            shading(ax2,'interp');colorbar('peer',ax2,'Location','Southoutside');  % Need to specify axes here.
            caxis(ax2,[-1 0.7]);
        else
            pcolor(ax2, ...
                       repmat(X_LOC:482,40,1), ...
                       squeeze(Z_REAL(X_LOC:482,Y_LOC,:))', ...
                       squeeze(VAR2_DEPTH(:,:,curr_date)'));
            shading(ax2,'interp');colorbar('peer',ax2,'Location','Southoutside');  % Need to specify axes here
            caxis(ax2,'auto')
        end
        % default to 200 m depth
        ax1.YLim =[-200 0];
        ax2.YLim = [-200 0];
        handles.edit_depthMax.String = ax1.YLim(1);
        handles.edit_depthMin.String = ax1.YLim(2);

        % Set color scale
        handles.edit_depth1ColorbarMin.String = round(ax1.CLim(1));       
        handles.edit_depth2ColorbarMin.String = round(ax2.CLim(1));
        if strncmp(VAR1_SELECTED,'Log10',5)
            handles.edit_depth1ColorbarMax.String = round(ax1.CLim(2),1);
        else
            handles.edit_depth1ColorbarMax.String = round(ax1.CLim(2));
        end
        if strncmp(VAR2_SELECTED,'Log10',5)
            handles.edit_depth2ColorbarMax.String = round(ax2.CLim(2),1);
        else
            handles.edit_depth2ColorbarMax.String = round(ax2.CLim(2));
        end        
        
        ax2.YTickLabels = '';  % Turn off for clarity
        
        colormap(ax1,jet);
        colormap(ax2,jet);
    else
        % 1 vertical plot
        ax = handles.axes_1vert;
        VAR1_SELECTED
        if strncmp(VAR1_SELECTED,'Log10',5)
            % Log10 Transform 
            pcolor(ax, ...
                       repmat(X_LOC:482,40,1), ...
                       squeeze(Z_REAL(X_LOC:482,Y_LOC,:))', ...
                       squeeze(log10(VAR1_DEPTH(:,:,curr_date)')));
            shading(ax,'interp');colorbar('peer',ax);
            caxis(ax,[-1 0.7]);
        else
            pcolor(ax, ...
                       repmat(X_LOC:482,40,1), ...
                       squeeze(Z_REAL(X_LOC:482,Y_LOC,:))', ...
                       squeeze(VAR1_DEPTH(:,:,curr_date)'));
            shading(ax,'interp');colorbar('peer',ax);
            caxis(ax,'auto');                        
        end
        % default to 200 m depth
        ax.YLim =[-200 0];
        handles.edit_depthMax.String = ax.YLim(1);
        handles.edit_depthMin.String = ax.YLim(2);

        % Set color scale
        handles.edit_depth1ColorbarMin.String = round(ax.CLim(1));        
        if strncmp(VAR1_SELECTED,'Log10',5)
            handles.edit_depth1ColorbarMax.String = round(ax.CLim(2),1);
        else
            handles.edit_depth1ColorbarMax.String = round(ax.CLim(2));
        end
        colormap(ax,jet);
    end
    

function edit_gridPointX_Callback(hObject, eventdata, handles)
% hObject    handle to edit_gridPointX (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_gridPointX as text
%        str2double(get(hObject,'String')) returns contents of edit_gridPointX as a double


% --- Executes during object creation, after setting all properties.
function edit_gridPointX_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_gridPointX (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


function edit_gridPointY_Callback(hObject, eventdata, handles)
% hObject    handle to edit_gridPointY (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_gridPointY as text
%        str2double(get(hObject,'String')) returns contents of edit_gridPointY as a double



% --- Executes during object creation, after setting all properties.
function edit_gridPointY_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_gridPointY (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


function edit_depthMin_Callback(hObject, eventdata, handles)
% hObject    handle to edit_depthMin (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_depthMin as text
%        str2double(get(hObject,'String')) returns contents of edit_depthMin as a double
    h_plot2button = handles.radiobutton_plot2;
    
    if h_plot2button(1).Value == 0 && h_plot2button(2).Value == 1
        % 2 vertical plots
        handles.axes_2vert1.YLim(2) = str2double(hObject.String);
        handles.axes_2vert2.YLim(2) = str2double(hObject.String);
    else
        % 1 vertical plot
        handles.axes_1vert.YLim(2) = str2double(hObject.String);
    end
    

% --- Executes during object creation, after setting all properties.
function edit_depthMin_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_depthMin (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
    hObject.String = 'Min';


function edit_depthMax_Callback(hObject, eventdata, handles)
% hObject    handle to edit_depthMin (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_depthMin as text
%        str2double(get(hObject,'String')) returns contents of edit_depthMin as a double
    h_plot2button = handles.radiobutton_plot2;
    
    if h_plot2button(1).Value == 0 && h_plot2button(2).Value == 1
        % 2 vertical plots
        handles.axes_2vert1.YLim(1) = str2double(hObject.String);
        handles.axes_2vert2.YLim(1) = str2double(hObject.String);
    else
        % 1 vertical plot
        handles.axes_1vert.YLim(1) = str2double(hObject.String);
    end

    
% --- Executes during object creation, after setting all properties.
function edit_depthMax_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_depthMin (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
    hObject.String='Max';

    
function edit_depth1ColorbarMin_Callback(hObject, eventdata, handles)
% hObject    handle to edit_depth1ColorbarMin (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_depth1ColorbarMin as text
%        str2double(get(hObject,'String')) returns contents of edit_depth1ColorbarMin as a double
    h_plot2button = handles.radiobutton_plot2;  % 2 vertical plots
    
    if h_plot2button(1).Value == 0 && h_plot2button(2).Value == 1
        % 2 horizontal plots
        handles.axes_2vert1.CLim(1) = str2double(hObject.String);
    else
        % 1 horizontal plot
        handles.axes_1vert.CLim(1) = str2double(hObject.String);
    end

% --- Executes during object creation, after setting all properties.
function edit_depth1ColorbarMin_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_depth1ColorbarMin (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit_depth1ColorbarMax_Callback(hObject, eventdata, handles)
% hObject    handle to edit_depth1ColorbarMax (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_depth1ColorbarMax as text
%        str2double(get(hObject,'String')) returns contents of edit_depth1ColorbarMax as a double
    h_plot2button = handles.radiobutton_plot2;  % 2 vertical plots
    
    if h_plot2button(1).Value == 0 && h_plot2button(2).Value == 1
        % 2 horizontal plots
        handles.axes_2vert1.CLim(2) = str2double(hObject.String);
    else
        % 1 horizontal plot
        handles.axes_1vert.CLim(2) = str2double(hObject.String);
    end
    

% --- Executes during object creation, after setting all properties.
function edit_depth1ColorbarMax_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_depth1ColorbarMax (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


function edit_depth2ColorbarMin_Callback(hObject, eventdata, handles)
% hObject    handle to edit_depth2ColorbarMin (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_depth2ColorbarMin as text
%        str2double(get(hObject,'String')) returns contents of edit_depth2ColorbarMin as a double
    handles.axes_2vert2.CLim(1) = str2double(hObject.String);


% --- Executes during object creation, after setting all properties.
function edit_depth2ColorbarMin_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_depth2ColorbarMin (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


function edit_depth2ColorbarMax_Callback(hObject, eventdata, handles)
% hObject    handle to edit_depth2ColorbarMax (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_depth2ColorbarMax as text
%        str2double(get(hObject,'String')) returns contents of edit_depth2ColorbarMax as a double
handles.axes_2vert2.CLim(2) = str2double(hObject.String);


% --- Executes during object creation, after setting all properties.
function edit_depth2ColorbarMax_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_depth2ColorbarMax (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes when selected object is changed in uibuttongroup_numPlots.
function uibuttongroup_numPlots_SelectionChangedFcn(hObject, eventdata, handles)
% hObject    handle to the selected object in uibuttongroup_numPlots 
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
    % Control visibility of plots and variable pop up menus based on number
    % of 2D-horizontal plots chosen.
    switch eventdata.NewValue.Tag
        case 'radiobutton_plot1'  % Plot 1 variable
            handles.popupmenu_var2.Visible = 'Off';
            handles.uipanel_caxis_plot2.Visible = 'Off';
            handles.uipanel_depth2Colorbar.Visible = 'Off';
            handles.axes_1horz.Visible = 'On';
            cla(handles.axes_1horz, 'reset');  % Refresh plot
            handles.axes_1vert.Visible = 'On';
            % reset has to be before clear, otherwise visibility is turned back on.
            cla(handles.axes_2horz1,'reset');  
            cla(handles.axes_2horz2,'reset');
            cla(handles.axes_2vert1,'reset');
            cla(handles.axes_2vert2,'reset');
            % Remove colorbars
            colorbar(handles.axes_2horz1,'hide');  
            colorbar(handles.axes_2horz2,'hide');
            colorbar(handles.axes_2vert1,'hide');
            colorbar(handles.axes_2vert2,'hide');            
            handles.axes_2horz1.Visible = 'Off';
            handles.axes_2horz2.Visible = 'Off';            
            handles.axes_2vert1.Visible = 'Off';
            handles.axes_2vert2.Visible ='Off';
             % Reset color axes
            handles.edit_caxis_plot1_min.String='Min';
            handles.edit_caxis_plot1_max.String='Max';
            disp('1 plot')
        case 'radiobutton_plot2'  % Plot 2 variables
            handles.popupmenu_var2.Visible = 'On';
            handles.uipanel_caxis_plot2.Visible = 'On';
            handles.axes_2horz1.Visible = 'On';
            handles.axes_2horz2.Visible = 'On';
            cla(handles.axes_1horz, 'reset');  % reset has to be before clear, otherwise visibility is turned back on.
            cla(handles.axes_1vert, 'reset'); 
            handles.axes_1horz.Visible = 'Off';
            handles.axes_2vert1.Visible = 'On';
            handles.axes_2vert2.Visible ='On';
            handles.axes_1vert.Visible = 'Off';
            handles.uipanel_depth2Colorbar.Visible = 'On';
            % Reset color axes
            handles.edit_caxis_plot1_min.String='Min';
            handles.edit_caxis_plot2_min.String='Min';
            handles.edit_caxis_plot1_max.String='Max';
            handles.edit_caxis_plot2_max.String='Max';
            disp('2 plots')
    end
    
    
% --- Executes during object creation, after setting all properties.    
function uibuttongroup_numPlots_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_depthMin (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
    
    
% --- Executes when selected object is changed in uibuttongroup_vertPlots.
function uibuttongroup_vertPlots_SelectionChangedFcn(hObject, eventdata, handles)
% hObject    handle to the selected object in uibuttongroup_vertPlots 
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

    % Control visibility of vertical plots and variable pop up menus based on
    % 2D-horizontal plots chosen.  Can turn off just the panel instead
    % of having to turn off individual plots, which allow the number of
    % vertical panels to also be controlled by the "Num Surf Plots"
    % button.
    switch eventdata.NewValue.Tag
        case 'radiobutton_vertPlotY'
            %handles.popupmenu_var2.Visible = 'On';
            handles.uipanel_vertical.Visible = 'On';
            handles.uipanel_paramVertical.Visible = 'On';
            handles.togglebutton_play.Visible = 'Off';
            disp('Vertical Plot')
        case 'radiobutton_vertPlotN' 
            handles.uipanel_vertical.Visible = 'Off';
            handles.uipanel_paramVertical.Visible = 'Off';
            handles.togglebutton_play.Visible = 'On';
    end
    
    
% --- Executes during object creation, after setting all properties.
function uibuttongroup_vertPlots_CreateFcn(hObject, eventdata, handles)
% hObject    handle to uibuttongroup_vertPlots (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called


% --- Executes during object creation, after setting all properties.
function radiobutton_vertPlotN_CreateFcn(hObject, eventdata, handles)
% hObject    handle to radiobutton_vertPlotN (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called
    hObject.Value = 1;  % Default to horizontal plot only.


% --- Executes during object creation, after setting all properties.    
function axes_1horz_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_depthMin (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes during object creation, after setting all properties.
function axes_2horz1_CreateFcn(hObject, eventdata, handles)
% hObject    handle to axes_2horz1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: place code in OpeningFcn to populate axes_2horz1
    hObject.Visible = 'Off';  % default to show only 1 plot.
    
    
    % --- Executes during object creation, after setting all properties.
function axes_2horz2_CreateFcn(hObject, eventdata, handles)
% hObject    handle to axes_2horz2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: place code in OpeningFcn to populate axes_2horz2
    hObject.Visible = 'Off';


% --- Executes during object creation, after setting all properties.
function axes_1vert_CreateFcn(hObject, eventdata, handles)
% hObject    handle to axes_1vert (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: place code in OpeningFcn to populate axes_1vert


% --- Executes during object creation, after setting all properties.
function axes_2vert1_CreateFcn(hObject, eventdata, handles)
% hObject    handle to axes_2vert1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: place code in OpeningFcn to populate axes_2vert1
hObject.Visible = 'Off';  % default to show only 1 plot.


% --- Executes during object creation, after setting all properties.
function axes_2vert2_CreateFcn(hObject, eventdata, handles)
% hObject    handle to axes_2vert2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: place code in OpeningFcn to populate axes_2vert2
hObject.Visible = 'Off';  % default to show only 1 plot.


% --- Executes during object creation, after setting all properties.
function uipanel_paramVertical_CreateFcn(hObject, eventdata, handles)
% hObject    handle to uipanel_paramVertical (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called
    hObject.Visible = 'Off';  % default to only show horizontal plot.


% --- Executes during object creation, after setting all properties.
function uipanel_vertical_CreateFcn(hObject, eventdata, handles)
% hObject    handle to uipanel_vertical (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called
hObject.Visible = 'Off';  % default to only show horizontal plot


% --- Executes during object creation, after setting all properties.    
function radiobutton_plot1_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_depthMin (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit_caxis_plot1_min_Callback(hObject, eventdata, handles)
% hObject    handle to edit_caxis_plot1_min (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_caxis_plot1_min as text
%        str2double(get(hObject,'String')) returns contents of edit_caxis_plot1_min as a double
    h_plot2button = handles.radiobutton_plot2;
    
    if h_plot2button(1).Value == 0 && h_plot2button(2).Value == 1
        % 2 horizontal plots
        handles.axes_2horz1.CLim(1) = str2double(hObject.String);
    else
        % 1 horizontal plot
        handles.axes_1horz.CLim(1) = str2double(hObject.String);
    end
    

% --- Executes during object creation, after setting all properties.
function edit_caxis_plot1_min_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_caxis_plot1_min (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


function edit_caxis_plot1_max_Callback(hObject, eventdata, handles)
% hObject    handle to edit_caxis_plot1_max (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_caxis_plot1_max as text
%        str2double(get(hObject,'String')) returns contents of edit_caxis_plot1_max as a double

    h_plot2button = handles.radiobutton_plot2;
    
    if h_plot2button(1).Value == 0 && h_plot2button(2).Value == 1
        % 2 horizontal plots
        handles.axes_2horz1.CLim(2) = str2double(hObject.String);
    else
        % 1 horizontal plot
        handles.axes_1horz.CLim(2) = str2double(hObject.String);
    end
    

% --- Executes during object creation, after setting all properties.
function edit_caxis_plot1_max_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_caxis_plot1_max (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


function edit_caxis_plot2_min_Callback(hObject, eventdata, handles)
% hObject    handle to edit_caxis_plot2_min (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_caxis_plot2_min as text
%        str2double(get(hObject,'String')) returns contents of edit_caxis_plot2_min as a double
    
    handles.axes_2horz2.CLim(1) = str2double(hObject.String);


% --- Executes during object creation, after setting all properties.
function edit_caxis_plot2_min_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_caxis_plot2_min (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


function edit_caxis_plot2_max_Callback(hObject, eventdata, handles)
% hObject    handle to edit_caxis_plot2_max (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_caxis_plot2_max as text
%        str2double(get(hObject,'String')) returns contents of edit_caxis_plot2_max as a double
     handles.axes_2horz2.CLim(2) = str2double(hObject.String);


% --- Executes during object creation, after setting all properties.
function edit_caxis_plot2_max_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_caxis_plot2_max (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes during object creation, after setting all properties.
function uipanel_caxis_plot2_CreateFcn(hObject, eventdata, handles)
% hObject    handle to uipanel_caxis_plot2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called
    hObject.Visible = 'Off';  % Default to 1 plot


%----------------------- Custom Function calls -----------------------%
function flist = get_flist_averages(exp_dir)
%
% Get filenames from selected experiment directory.
% Input: exp_dir: cell array of the directory to read files from
% Output: flist: cell array of all NetCDF average files in exp_dir.
%
    header = '/home/morrison/data2/zach/';
    d = dir(strcat(header,string(exp_dir),'/*avg*'));
    nfiles = length({d.name});
    flist={nfiles};  % preallocated for speed
    
    if nfiles == 0
        disp('Empty Directory')
        return
    elseif nfiles == 1
        % Ensure file is represented as a cell.  
        flist = cellstr(strcat(d.folder,'/',d.name));
    else
        for i = 1:nfiles
            flist{i} = strcat(d(i).folder,'/',d(i).name);
        end
    end

    return
    
function vartime = combine_var_times(flist, varname, depth_switch) 
%
% Generate time series of selected variable.
% Input: flist (cell array): files from which to read 'varname'.
%           varname (string): variable of interest to plot.
%           depth_switch(bool): toggle if start_3D and count_3D handle
%                                           depth data.
% Output: vartime: time series of variable 'varname'.
%
    global X_LOC Y_LOC
    
    %Initialize start and end points for 3D data slices.
    if depth_switch
        % Read depth data for variable at X and Y locations specified by
        % user.  Read X-data to boundary for cross-shelf view.
        start_3D = [X_LOC,Y_LOC,1,1];
        count_3D = [Inf,1,Inf,Inf];
    else
        % Read surface data
        start_3D = [1,1,40,1];
        count_3D = [Inf,Inf,Inf,Inf];
    end
    
    % Read in first file to determine variable rank
    if strcmp(varname,'Total Chl') || strcmp(varname,'Log10 Total Chl')         
        % Combine Chlorophyll terms
        vartime = squeeze(ncread(flist{1},'Chlorophyll',start_3D,count_3D) ...
                      + ncread(flist{1},'ChlD',start_3D,count_3D));
    elseif strcmp(varname,'Log10 Total Chl 1.59')
        % constant Chl:N ratio
        vartime = 1.59*squeeze(ncread(flist{1},'Phytoplankton',start_3D,count_3D) ...
                     + ncread(flist{1},'Diatoms',start_3D,count_3D));
    elseif strcmp(varname,'Total Phy') || strcmp(varname,'Log10 Total Phy')
        vartime = squeeze(ncread(flist{1},'Phytoplankton',start_3D,count_3D) ...
                      + ncread(flist{1},'Diatoms',start_3D,count_3D));
    elseif strcmp(varname,'Log10 ChlS')
        vartime = squeeze(ncread(flist{1},'Chlorophyll',start_3D,count_3D));
    elseif strcmp(varname,'Log10 ChlS 1.59')
        % Constant Chl:N for SPHY
        vartime = 1.59*squeeze(ncread(flist{1},'Phytoplankton',start_3D,count_3D));
    elseif strcmp(varname,'Log10 ChlD')
        vartime = squeeze(ncread(flist{1},'ChlD',start_3D,count_3D));
    elseif strcmp(varname,'Log10 ChlD 1.59')
        % Constant Chl:N for LPHY
        vartime = 1.59*squeeze(ncread(flist{1},'Diatoms',start_3D,count_3D));
    elseif strncmp(varname,'Chl:N',5)
        % Chl:N ratios
        switch varname
            case 'Chl:N Tot'
                vartime = squeeze(ncread(flist{1},'Chlorophyll',start_3D,count_3D) ...
                             + ncread(flist{1},'ChlD',start_3D,count_3D)) ./ ...
                                 squeeze(ncread(flist{1},'Phytoplankton',start_3D,count_3D) ...
                             + ncread(flist{1},'Diatoms',start_3D,count_3D));
            case 'Chl:N S'
                 vartime = squeeze(ncread(flist{1},'Chlorophyll',start_3D,count_3D)) ./ ...
                                  squeeze(ncread(flist{1},'Phytoplankton',start_3D,count_3D));
            case 'Chl:N L'
                vartime = squeeze(ncread(flist{1},'ChlD',start_3D,count_3D)) ./ ...
                                  squeeze(ncread(flist{1},'Diatoms',start_3D,count_3D));
        end
    elseif strcmp(varname, 'ocean_time')
        vartime = ncread(flist{1}, varname);
    else
        % Normal variable
        vartime = squeeze(ncread(flist{1}, varname,start_3D,count_3D));
    end
    
    if ~exist('h','var')
        global h
        h = ncread(flist{1},'h');
    end
    ndims = size(size(vartime));
    
    % Concatenate tracer along appropriate time dimension.
    switch ndims(2)
        case 2  % scalar variable, e.g. ocean_time, time is second dimension
            if isempty(flist)
                disp('Empty Directory')
                return
            elseif length(flist) == 1
                % already read in only available file.
                %vartime
                return  
            else
                for i=2:length(flist)
                vartime = cat(1,vartime,ncread(flist{i},varname));
                end
            end
       case 3 % 2D or 3D variable, doesn't matter since we load only a single vertical level (for speed)
            if isempty(flist)
                disp('Empty Directory')
                return
            elseif length(flist) == 1
                % already read in only available file.
                %vartime
                return
            else
                for i=2:length(flist)
                    if strcmp(varname,'Total Chl') || strcmp(varname,'Log10 Total Chl')
                    % Combine Chlorophyll terms
                        vartime = cat(3, vartime, ...
                                               squeeze(ncread(flist{i},'Chlorophyll',start_3D,count_3D) ...
                                      +       ncread(flist{i},'ChlD',start_3D,count_3D)));
                    elseif strcmp(varname,'Log10 Total Chl 1.59')
                        % constant Chl:N ratio
                        vartime = cat(3, vartime, ...
                                               1.59*squeeze(ncread(flist{i},'Phytoplankton',start_3D,count_3D) ...
                                     +        ncread(flist{i},'Diatoms',start_3D,count_3D)));
                    elseif strcmp(varname,'Total Phy') || strcmp(varname,'Log10 Total Phy')
                        % Combine phyto and diatom terms
                        vartime = cat(3, vartime, ...
                                               squeeze(ncread(flist{i},'Phytoplankton',start_3D,count_3D) ...
                                      +      ncread(flist{i},'Diatoms',start_3D,count_3D)));
                    elseif strcmp(varname,'Log10 ChlS')
                        vartime = cat(3,vartime,squeeze(ncread(flist{i},'Chlorophyll',start_3D,count_3D)));
                    elseif strcmp(varname,'Log10 ChlS 1.59')
                        % constant Chl:N ratio for SPHY
                        vartime = cat(3, vartime, ...
                                               1.59*squeeze(ncread(flist{i},'Phytoplankton',start_3D,count_3D)));
                    elseif strcmp(varname,'Log10 ChlD')
                        vartime = cat(3,vartime,squeeze(ncread(flist{i},'ChlD',start_3D,count_3D)));
                    elseif strcmp(varname,'Log10 ChlD 1.59')
                        % constant Chl:N ratio for LPHY
                        vartime = cat(3, vartime, ...
                                               1.59*squeeze(ncread(flist{i},'Diatoms',start_3D,count_3D)));
                    elseif strncmp(varname,'Chl:N',5)
                        % Chl:N ratios
                        switch varname
                            case 'Chl:N Tot'
                                vartime = cat(3,vartime, ...
                                                       squeeze(ncread(flist{i},'Chlorophyll',start_3D,count_3D) ...
                                             +        ncread(flist{i},'ChlD',start_3D,count_3D)) ./ ...
                                                        squeeze(ncread(flist{i},'Phytoplankton',start_3D,count_3D) ...
                                             +        ncread(flist{i},'Diatoms',start_3D,count_3D)));
                            case 'Chl:N S'
                                 vartime = cat(3,vartime, ...
                                                        squeeze(ncread(flist{i},'Chlorophyll',start_3D,count_3D)) ./ ...
                                                        squeeze(ncread(flist{i},'Phytoplankton',start_3D,count_3D)));
                            case 'Chl:N L'
                                vartime = cat(3,vartime, ...
                                                       squeeze(ncread(flist{i},'ChlD',start_3D,count_3D)) ./ ...
                                                       squeeze(ncread(flist{i},'Diatoms',start_3D,count_3D)));
                        end                                           
                    else
                        vartime = cat(3,vartime,squeeze(ncread(flist{i},varname,start_3D,count_3D)));
                    end
                end
            end           
    end
    
    return 

function grid_container = get_grid_vars()
%
% Read in grid variables relevant for plotting purposes.
% Input: None.
% Output: lat: latitude at rho-points
%              lon: longitude at rho-points
%              N: Number of vertical levels
% 

    global h
    % Establish grid file to read from.
    f1 = '/home/server/homes/pi/zwallace/';
    f2 = 'ROMS/Project_Patagonia/bulk';
    f3 = 'roms_grd_rivers.nc';
    grdfile = fullfile(f1,f2,f3);
    
    if ~exist('h','var')
        h=ncread(grdfile,'h');
    end
    
    % Build cell array to loop over and grab variables.  Will help if
    % wanting to quickly add new variables to read in the future.
    % Make sure 'N' is the last key, since it's not being read from the
    % grid file currently -- Z.W. 1.17.19 
    key_set = {'lat_rho',...
                         'lon_rho',...
                         'N'};
    value_set = {length(key_set)};
    
    % Open grid file for reading and obtain grid variables of interest.
    ncid = netcdf.open(grdfile);
    
    % s_rho is not in grid file, set to 40 for now
    for i = 1:length(key_set)-1
        varid = netcdf.inqVarID(ncid, key_set{i});
        value_set{i} = netcdf.getVar(ncid,varid,'double');
    end
    value_set{3} = 40;  % Append s_rho at end
    
    netcdf.close(ncid);
    
    % Build container with read grid variables.
     grid_container = containers.Map(key_set,value_set);
     
    return 


% ---------- FUNCTION HANDLES DEFINED B.C. ACCIDENTALLY DELETED ---------%
% --- Executes during object creation, after setting all properties.
function uipanel_lon_CreateFcn(hObject, eventdata, handles)
% hObject    handle to uipanel_lon (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% --- Executes during object creation, after setting all properties.
function edit_latmin_CreateFcn(hObject, eventdata, handles)
% hObject    handle to uipanel_lon (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% --- Executes during object creation, after setting all properties.
function edit_latmax_CreateFcn(hObject, eventdata, handles)
% hObject    handle to uipanel_lon (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% --- Executes during object creation, after setting all properties.
function edit_lonmin_CreateFcn(hObject, eventdata, handles)
% hObject    handle to uipanel_lon (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% --- Executes during object creation, after setting all properties.
function edit_lonmax_CreateFcn(hObject, eventdata, handles)
% hObject    handle to uipanel_lon (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% --- Executes during object creation, after setting all properties.
function figure1_CreateFcn(hObject, eventdata, handles)
% hObject    handle to uipanel_lon (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% --- Executes during object deletion
function uipanel_lon_DeleteFcn(hObject, eventdata, handles)
% hObject    handle to uipanel_lon (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

%---------------------------------------------------------------------------------------------%
