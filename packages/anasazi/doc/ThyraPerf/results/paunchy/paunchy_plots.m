Al = [0,0,0,0,0,0];
Ah = [850,75,950,310,140,625];
L = ['1e2';'1e3';'1e4';'1e5';'1e6'];
files = ...
 { ...
   { 'LOBPCG-S-LA.dat', 'LOBPCG-2-LA.dat',  'LOBPCG-4-LA.dat' } ...
   { 'LOBPCG-S-OP.dat', 'LOBPCG-2-OP.dat',  'LOBPCG-4-OP.dat' } ...
   { 'LOBPCG-S-TO.dat', 'LOBPCG-2-TO.dat',  'LOBPCG-4-TO.dat' } ...
   { 'BD-S-LA.dat',     'BD-2-LA.dat'    ,  'BD-4-LA.dat'     } ...
   { 'BD-S-OP.dat',     'BD-2-OP.dat'    ,  'BD-4-OP.dat'     } ...
   { 'BD-S-TO.dat',     'BD-2-TO.dat'    ,  'BD-4-TO.dat'     } ...
 };

machname = 'Paunchy';
codename = {'BD','BD','BD','LOBPCG','LOBPCG','LOBPCG'};
meas  = {'LinAlg', 'ApplyOp', 'Total','LinAlg', 'ApplyOp', 'Total'};
procs = {'S', '2', '4'};

nmeas = length(meas);
nprocs = length(procs);

for i=1:nmeas,
    for j=1:nprocs,
        fn = files{i}{j};
        M{i}{j} = load(fn);
    end
end

npts = size(M{1}{1},1);

p1 = linspace(.2,.8,npts);
for i=1:nmeas,
    figure;
    hold on;
    for j=1:nprocs,
        h1 = plot(j-1+p1,M{i}{j}(:,1),'*-r');
        h2 = plot(j-1+p1,M{i}{j}(:,2),'o--g');
        h3 = plot(j-1+p1,M{i}{j}(:,3),'+:b');
    end
    hold off;
    a = axis;
    axis([0,nprocs,Al(i),Ah(i)]);
    set(gca,'XGrid','off');
    set(gca,'XMinorGrid','off');
    set(gca,'YGrid','on');
    set(gca,'YMinorGrid','off');
    set(gca,'FontSize',12);
    set(gca,'XTickLabel',procs);
    set(gca,'XTickMode','manual');
    set(gca,'XTick',.5+[0:nprocs-1]);
    xlabel('# processors');
    ylabel('Runtime (s)');
    title(sprintf('%s - %s - %s',machname,codename{i},meas{i}));
    legend([h1 h2 h3],{'Epetra','Thyra','Thyra-Epetra'},'Location','NorthEast');
    p = get(gcf,'Position');
    p(4) = p(4)*.75;
    set(gcf,'Position',p);

    saveas(gcf,sprintf('%s-%s-%s_ln.jpg',machname,codename{i},meas{i}),'jpeg');
    saveas(gcf,sprintf('%s-%s-%s_ln.eps',machname,codename{i},meas{i}),'epsc');
    refresh(gcf);
end
