function plot_sim(sim)
fig = gcf;
if isobject(fig), fig = fig.Number; end

Trait_names = {sim.Traits(:).name};
Trait_dims = [sim.Traits(:).dims];

% Remove constants:
lociPD = [sim.Traits(:).loci_per_dim];
Trait_names(lociPD==0) = [];
Trait_dims(lociPD==0) = [];

switch lower(sim.Space.model)
    case 'none'
    case {'discrete','continuous'}
        Trait_names = {Trait_names{:},'pos'};
        Trait_dims = [Trait_dims sim.Space.dimensions];
    otherwise
        error(['Unknown Space model : ' sim.Space.model])
end
nT = length(Trait_names);
        
% find plot ranges all traits:
bins = cell(1,nT);
if strcmpi(sim.Genetics.model,'diallelic')
    maxx = zeros(1,nT);
    for tr = 1:nT
        x = [sim.samples.(Trait_names{tr})];
        maxx(tr) = max(abs(x(:)));
    end
    
    for si=1:length(sim.samples)
        s = sim.samples(si);
        for tr = 1:nT
            x = s.(Trait_names{tr});
            if maxx(tr) > 0
                x = round(x/maxx(tr),3)*maxx(tr);
            end
            bins{tr} = unique([bins{tr} x(:)']);
        end
    end
    for tr=1:nT
        if length(bins{tr})>1
            bins{tr} = min(bins{tr}):min(diff(bins{tr})):max(bins{tr});
        else
            bins{tr} = bins{tr}+[-1 0 1];
        end
    end
else
    Smaxv = -inf + zeros(1,length(Trait_names));
    Sminv = inf + zeros(1,length(Trait_names));
    for i=1:length(sim.samples)
        s = sim.samples(i);
        for t = 1:length(Trait_names)
            x = s.traits.(Trait_names{t});
            m = max(x(:));
            Smaxv(t) = max(Smaxv(t),m);
            m = min(x(:));
            Sminv(t) = min(Sminv(t),m);
        end
    end
    for t=1:nT
        if Smaxv(t)>Sminv(t)
            bins{t} = linspace(Sminv(t),Smaxv(t),50);
        else
            bins{t} = linspace(Sminv(t)-1,Smaxv(t)+1,3);
        end
    end
end
% plot all traits:
figure(fig), clf
tt = [sim.samples.gen];
nn = [sim.samples.size];
%sp = 1;
rows = max(Trait_dims);
for tr=1:nT
    for d=1:Trait_dims(tr)
        positions = zeros(length(bins{tr}),sim.sample_count);
        for si = 1:sim.sample_count
            s = sim.samples(si);
            if s.size > 0
                xx = s.(Trait_names{tr});
                xx = xx(d,:);
                if 0 %isfield(s,'pos')
                    pp = s.pos;
                    % mean position in each bin:
                    for p = 0:max(pp)
                        Np = hist(xx(pp==p),bins{tr})';
                        positions(:,si) = positions(:,si) + p*Np;
                    end
                    N = hist(xx,bins{tr})';
                    positions(N>0,si) = positions(N>0,si)./N(N>0);
                    positions(N==0,si) = NaN;
                else
                    positions(:,si) = hist(xx,bins{tr})';
                    positions(positions(:,si)==0,si) = NaN;
                end
            end
        end
        subplot(rows,nT,(d-1)*nT+tr);
        pcolorC(tt,bins{tr},positions)
        shading flat
        xlabel('generation')
        if Trait_dims(tr)>1
            ylabel([Trait_names{tr} '_' num2str(d)])
        else
            ylabel(Trait_names{tr})
        end
    end
end

figtitle(sim.name)

% figure(fig+1), clf
% %[mean_Achoose,mean_Bchoose,std_colour,nnsp] = plot_RI_trends(sim,false,0.1);
% subplot(2,1,1)
% if tt(end) < sim.t_max % extinction?
%     nn = [nn 0];
%     tt = [tt tt(end)+sim.sample_interval];
% end
% plot(tt,nn,'b-')
% xlabel('generations')
% ylabel('total population size (blue)')
% raxis
% plot(tt,nnsp,'r-');
% ylabel('#species (red)')
% ylim([0 max(max(nnsp),20)])

% subplot(2,1,2)
% plot(tt, mean_Achoose,'b', tt,mean_Bchoose,'g')
% xlabel('generations')
% ylabel('mean(Ac)(blue), mean(Bc)(red)')
% ylim([0 10])
% raxis
% plot(tt, std_colour,'r');
% ylabel('SD(colour)(red)')
% ylim([0 20])

figtitle(sim.name)
drawnow

return
figure(fig+2), clf
i1 = 0;
nD = sum(Trait_dims);
sa = sim.samples(end);
for tr1=1:nT
    for d1=1:Trait_dims(tr1)
        i1 = i1+1;
        i2 = 0;
        for tr2=1:nT
            for d2=1:Trait_dims(tr2)
                i2 = i2+1;
                subplot(nD,nD,(i1-1)*nD + i2)
                if i1==i2
                    histogram(sa.(Trait_names{tr1})(d1,:), 50, 'DisplayStyle', 'Stairs')
                    xlabel([Trait_names{tr2} '_' num2str(d2)])
                else
                    H = histogram2(sa.(Trait_names{tr2})(d2,:), sa.(Trait_names{tr1})(d1,:), 20, 'DisplayStyle','tile');
                    H.EdgeColor = 'none';
                    grid off
                    corr = corrcoef(sa.(Trait_names{tr2})(d2,:), sa.(Trait_names{tr1})(d1,:));
                    title(sprintf('\\rho=%0.3f', corr(1,2)))
                    if i1==nD
                        xlabel([Trait_names{tr2} '_' num2str(d2)])
                    end
                    if i2==1
                        ylabel([Trait_names{tr1} '_' num2str(d1)])
                    end
                end
            end
        end
    end
end
figtitle(['Last sample distributions, ' sim.name])

drawnow

function pcolorC(X,Y,C)
%centered pcolor
if min(size(X))>1 % matrix input?
    X=X(1,:);
end
X = X(:)';
% Find midpoints between X-values. They will be gridlines.
if length(X)>1
    xmid = (X(1:end-1) + X(2:end))/2;
    % First point to the left:
    x0 = X(1) - (xmid(1) - X(1));
    % And last to the right:
    xend = X(end) + (X(end)-xmid(end));
    xx = [x0 xmid xend];
else
    xx = [X-0.5 X+0.5];
end
    
if min(size(Y))>1 % matrix input?
    Y=Y(:,1);
end

Y = Y(:)';
% Find midpoints between Y-values. They will be gridlines.
if length(Y)>1
ymid = (Y(1:end-1) + Y(2:end))/2;
% First point to the left:
y0 = Y(1) - (ymid(1) - Y(1));
% And last to the right:
yend = Y(end) + (Y(end)-ymid(end));
yy = [y0 ymid yend];
else
    yy = [Y-0.5 Y+0.5];
end
% Pad with NaNs:
C(end+1,:) = NaN;
C(:,end+1) = NaN;
pcolor(xx,yy,C);

function figtitle(tit)
ca = get(gcf,'currentaxes');
ha = axes('position',[0.5 1 0.1 0.1]); % [left, bottom, width, height]
set(ha,'visible','off','userdata','pagetitle')
set(ha,'units','characters')
pos = get(ha,'position');
set(ha,'position',[pos(1) pos(2)-1 pos(3) 1])
text(0,0,tit,'horizontalAlignment','Center','fontsize',12);
set(ha,'units','normalized')
set(gcf,'currentaxes',ca)

