clear all
start_iter = 250000;
end_iter   = 350000;
step_iter  = 2000;
nz=384;
ns=floor((end_iter-start_iter)/step_iter)+1;

s=zeros(nz,13,ns);

counter=0;
for iter = start_iter:step_iter:end_iter
    filename = sprintf('stat_%08d.dat', iter);
    counter=counter+1;
    s(:,:,counter) = load(filename);
end

for i=1:ns
    means(:,:)=mean(s,3);
end

mean180=importdata('chan180/chan180.means');
rey180=importdata('chan180/chan180.reystress');

figure(1)
clf
% mean velocity profile
subplot(2,2,1)
plot(s(:,1),means(:,2),'LineWidth',3.,'DisplayName','TCF36')
hold on
plot(mean180(:,1),mean180(:,3),'-r','LineWidth',3.,'DisplayName','Moser');
legend show
hold off

% RMS u (streamwise)
subplot(2,2,2)
plot(s(:,1),means(:,5),'LineWidth',3.,'DisplayName','TCF36')
hold on
plot(mean180(:,1),sqrt(rey180(:,3)),'-r','LineWidth',3.,'DisplayName','Moser');
ylim([0 3]);
legend show
hold off

% RMS v (spanwise)
subplot(2,2,3)
plot(s(:,1),means(:,6),'LineWidth',3.,'DisplayName','TCF36')
hold on
plot(mean180(:,1),sqrt(rey180(:,5)),'-r','LineWidth',3.,'DisplayName','Moser');
ylim([0 3]);
legend show
hold off

% RMS w (wall-normal)
subplot(2,2,4)
plot(s(:,1),means(:,7),'LineWidth',3.,'DisplayName','TCF36')
hold on
plot(mean180(:,1),sqrt(rey180(:,4)),'-r','LineWidth',3.,'DisplayName','Moser');
ylim([0 3]);
legend show
hold off