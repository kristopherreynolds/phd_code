classdef hermite<handle
    %hermite functions class
    %Author: Kristopher Reynolds
    properties
        P %parameter structure
        fbar %where the discrete data lives
        Hfuns %hermite functions matrix
        fhats %hermite coefficients
        fherm %bandlimited hermite expansion
        nlse %normalized least-squares error (%)
        min_sval %minimum singular value
    end
    
    methods
        function obj = hermite(P,fbar)
            obj.P = P;
            obj.fbar = fbar;
        end
        
        function fit_series(obj)
            if ~isfield(obj.P,'estimation_method')
                obj.P.estimation_method = 'svd';
            end
            
            if numel(size(obj.fbar))==2 && min(size(obj.fbar))==1 %verify fbar is 1d
                N = obj.P.N;
                a = obj.P.a;
                npts = numel(obj.fbar);
                range = obj.P.range;
                x = linspace(-range,range,npts);
                obj.P.x = x; %store x in parameters structure
                obj.Hfuns = herm_funs(N,a*x);
                obj.fhats = zeros(N+1,1);
                if strcmp(obj.P.estimation_method,'svd')
                    %estimate 1d coefficients with SVD
                    bool1 = npts > N + 1;
                    if ~bool1
                        disp('SVD Invalid')
                        return
                    end
                    [U,S,V] = svd(obj.Hfuns);
                    S1 = S(1:N+1,:);
                    s = diag(S1);
                    obj.min_sval = min(s); %minimum singular value
                    U1 = U(:,1:N+1);
                    Hinv = V/S1*U1';
                    obj.fhats = Hinv*obj.fbar;
                    obj.fherm = obj.Hfuns*obj.fhats;
                    obj.nlse = 100*sum((obj.fbar-obj.fherm).^2)/sum((obj.fbar).^2);
                elseif strcmp(obj.P.estimation_method,'fft')
                    %estimate 1d coefficients with FFT fit
                    %COMING SOON
                end
            elseif numel(size(obj.fbar))==2 && min(size(obj.fbar))~=1 %verify fbar is 2d
                %note: assume data is square
                N = obj.P.N;
                a = obj.P.a;
                npts = size(obj.fbar,1);
                ii = 1 : npts;
                range = obj.P.range;
                x = -range + 2*range*(ii-1)/(npts-1); clear ii
                obj.P.x = x; %store x in parameters structure
                obj.Hfuns = herm_funs(N,a*x);
                if strcmp(obj.P.estimation_method,'svd')
                    %estimate 2d coefficients with SVD
                    bool1 = npts > N+1;
                    if ~bool1
                        disp('SVD Invalid')
                        return
                    end
                    [U,S,V] = svd(obj.Hfuns);
                    S1 = S(1:N+1,:);
                    s = diag(S1);
                    obj.min_sval = min(s); %minimum singular value
                    U1 = U(:,1:N+1);
                    Hinv = V/S1*U1';
                    obj.fhats = Hinv*obj.fbar*Hinv';
                    if isfield(obj.P,'special_limits')
                        if obj.P.special_limits
                            %we are going to enforce the
                            %Bandlimit-Invariance-to-Rotation Limits
                            obj.fherm = obj.fbar*0;
                            fhats2 = obj.fhats*0;
                            for m = 0 : N
                                for n = 0 : N-m
                                    fhats2(m+1,n+1) = obj.fhats(m+1,n+1);
                                end
                            end
                            obj.fhats = fhats2;
                        end
                        
                    end
                    obj.fherm = obj.Hfuns*obj.fhats*obj.Hfuns';
                    obj.nlse = 100*sum(sum((obj.fbar-obj.fherm).^2))/sum(sum((obj.fbar).^2));
                else
                    %estimate 2d coefficients with FFT fit
                    %COMING SOON
                end
            elseif numel(size(obj.fbar))==3
                %3D Hermite Series Fit
                %COMING SOON
                
            end
        end
        
        function plot(obj)
            if numel(size(obj.fbar))==2 && min(size(obj.fbar))==1 %verify fbar is 1d
                if ~isempty(obj.fherm)
                    figure
                    plot(obj.P.x,obj.fbar,obj.P.x,obj.fherm,'r--','linewidth',2)
                    xlabel('x')
                    ylabel('Data')
                    title(['NLSE = ',num2str(obj.nlse)])
                    grid on
                    legend('Data','Hermite Series')
                else
                    disp('Hermite Expansion not Defined')
                end
            elseif numel(size(obj.fbar))==2 && min(size(obj.fbar))~=1 %verify fbar is 2d
                if ~isempty(obj.fherm)
                    figure
                    subplot(1,2,1)
                    imagesc(obj.P.x,obj.P.x,obj.fbar)
                    %axis xy
                    colorbar
                    if isfield(obj.P,'clims')
                        caxis (obj.P.clims)
                    end
                    if isfield(obj.P,'cmap')
                        colormap(obj.P.cmap)
                    end
                    axis equal
                    xlabel('x')
                    ylabel('y')
                    subplot(1,2,2)
                    imagesc(obj.P.x,obj.P.x,abs(obj.fbar-obj.fherm))
                    %axis xy
                    colorbar
                    axis equal
                    xlabel('x')
                    ylabel('y')
                else
                    disp('Hermite Expansion not Defined')
                end
            end
        end
        
        
        
    end
    
end



















