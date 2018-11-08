datasetnames={'DBCG','Osl2','Micma'};
%% import DBCG
DBCG=importdata('C:\Users\shaybe\Dropbox\Thesis-PHd\QuantileNormalization\rawdata\DBCG_miRNA_raw_154samples_simple_with_headers_removed_missing.csv');
DBCGtbl=array2table(DBCG.data,'VariableNames',cellfun(@(x) horzcat('x',replace(replace(x,'.','_'),'-','_')),DBCG.textdata(1,3:end),'uniformoutput',false), 'RowNames',DBCG.textdata(2:end,2));

%% import OSL2
Osl2=importdata('C:\Users\shaybe\Dropbox\Thesis-PHd\QuantileNormalization\rawdata\OSL2_miRNA_raw_425_simple_with_headers_removed_missing.csv');
Osl2tbl=array2table(Osl2.data,'VariableNames',cellfun(@(x) horzcat('x',replace(replace(x,'.','_'),'-','_')),Osl2.textdata(1,3:end),'uniformoutput',false), 'RowNames',Osl2.textdata(2:end,2));

%% import Micma
Micma=importdata('C:\Users\shaybe\Dropbox\Thesis-PHd\QuantileNormalization\rawdata\Micma_miRNA_raw_200samples_in replicates_simple_with_headers_removed_missing.csv');
Micmatbl=array2table(Micma.data,'VariableNames',cellfun(@(x) horzcat('x',replace(replace(x,'.','_'),'-','_')),Micma.textdata(1,3:end),'uniformoutput',false), 'RowNames',Micma.textdata(2:end,2));
MicmaGroupings=importdata('C:\Users\shaybe\Dropbox\Thesis-PHd\QuantileNormalization\quantile_normalization_variant\main\reference_files\Experiment grouping file.csv');
MicmaClinicalLabels=importdata('C:\Users\shaybe\Dropbox\Thesis-PHd\QuantileNormalization\quantile_normalization_variant\main\cohorts_clinical_info\Micma_samples_clinical_info.csv');

%% Join data
numsmpls=[size(DBCG.data,2) size(Osl2.data,2) size(Micma.data,2)];
[joined1,ia1,ib1]=outerjoin(DBCGtbl,Osl2tbl,'Keys','Row');
[Mtbl,ia2,ib2]=outerjoin(joined1,Micmatbl,'Keys','Row');
M=table2array(Mtbl);

%% Normalize data
Mnorm = quantilenorm(M);
Mnorm2 = quantilenorm(M,'median',true);

DBCGNorm=Mnorm(:,1:numsmpls(1));
Osl2Norm=Mnorm(:,numsmpls(1)+1:numsmpls(1)+numsmpls(2));
MicmaNorm=Mnorm(:,end-numsmpls(3)+1:end);

%% Plot tsne impact
prplx=20;
goodrows=not(any(isnan(M),2));
[preY,preloss]=tsne(log(M(goodrows,:))','algorithm','exact','perplexity',prplx);
goodrows=not(any(isnan(Mnorm),2));
[postY,postloss]=tsne(log(Mnorm(goodrows,:))','algorithm','exact','perplexity',prplx);
goodrows=not(any(isnan(Mnorm2),2));
[postY2,postloss2]=tsne(log(Mnorm2(goodrows,:))','algorithm','exact','perplexity',prplx);

datasourceids=arrayfun(@(x,id) repmat(id,1,x), numsmpls, 1:length(numsmpls),'uniformoutput',false);
figure('color','w');
subplot(1,3,1); 
gscatter(preY(:,1),preY(:,2),datasetnames([datasourceids{:}])');
title(horzcat('Pre normalization $\mathcal{L}$=', num2str(preloss)),'Interpreter','latex')
subplot(1,3,2);
gscatter(postY(:,1),postY(:,2),datasetnames([datasourceids{:}])');
title(horzcat('Quantile normalization $\mathcal{L}$=', num2str(postloss)),'Interpreter','latex');
subplot(1,3,3);
gscatter(postY2(:,1),postY2(:,2),datasetnames([datasourceids{:}])');
title(horzcat('Quantile normalization (median) $\mathcal{L}$=', num2str(postloss2)),'Interpreter','latex');
suptitle(horzcat('tSNE (perplexity=',num2str(prplx),')'));

%% Clustering visualization
figure('color','w');
clustergram

%% Micma replicates comparison
coefdiff=nan(max(MicmaGroupings.data),2);
for gid=1:max(MicmaGroupings.data)
    if (sum(MicmaGroupings.data==gid) ~= 2)
        display(horzcat('replicate ', num2str(gid),' is of size ',num2str(sum(MicmaGroupings.data==gid))));
        continue;
    end
    [~,currSet]=intersect(Micma.textdata(1,3:end), MicmaGroupings.textdata([false;MicmaGroupings.data==gid],1));
    f=figure('color','w'); 
    scatter(Micma.data(:,currSet(1)),Micma.data(:,currSet(2)),'b'); hold on;
    scatter(MicmaNorm(:,currSet(1)),MicmaNorm(:,currSet(2)),'r');
    maxval=max(max([Micma.data(:,currSet); MicmaNorm(:,currSet)]));
    plot([1 maxval], [1 maxval],'--k')
    set(gca,'xscale','log','yscale','log');
    [op_rho,op_pval]=corr(Micma.data(:,currSet(1)),Micma.data(:,currSet(2)),'type','pearson','rows','complete');
    [os_rho,os_pval]=corr(Micma.data(:,currSet(1)),Micma.data(:,currSet(2)),'type','spearman','rows','complete');
    [np_rho,np_pval]=corr(MicmaNorm(:,currSet(1)),MicmaNorm(:,currSet(2)),'type','pearson','rows','complete');
    [ns_rho,ns_pval]=corr(MicmaNorm(:,currSet(1)),MicmaNorm(:,currSet(2)),'type','spearman','rows','complete');
    coefdiff(gid,:)=[np_rho-op_rho,ns_rho-os_rho];
    legend(horzcat('Original (pearson=',num2str(op_rho),' ,spearman=',num2str(os_rho),')'),...
        horzcat('QuantileNormalized (pearson=',num2str(np_rho),' ,spearman=',num2str(ns_rho),')'),'location','best')
    title(MicmaGroupings.textdata{1+find(MicmaGroupings.data==gid,1),2});
    [~,~,~]=mkdir(horzcat(pwd,'\Micma\ReplicatesComparison\'));
    saveas(f,horzcat(pwd,'\Micma\ReplicatesComparison\',MicmaGroupings.textdata{1+find(MicmaGroupings.data==gid,1),2},'.png'))
    close(f);
end 
%
f=figure('color','w');
hist(coefdiff(:,1),40); xlabel('pearson(replicates normalized)-pearson(replicates raw)');
[~,p]=ttest(coefdiff(:,1),0,'tail','left');
title({'Impact of normalization on replicate correlation'; horzcat('t-test left tail: p=',num2str(p),')')})


%% Micma embeddings
prplx=20;
goodrows=not(any(isnan(Micma.data),2));
[preY,preloss]=tsne(Micma.data(goodrows,:)','algorithm','exact','perplexity',prplx);
goodrows=not(any(isnan(MicmaNorm),2));
[postY,postloss]=tsne(MicmaNorm(goodrows,:)','algorithm','exact','perplexity',prplx);
%
figure('color','w');
[~,subset,ord]=intersect(Micma.textdata(1,3:end), MicmaGroupings.textdata(2:end,1));
cmap=colormap(lines);
subplot(2,5,1); 
gscatter(preY(:,1),preY(:,2),MicmaGroupings.data(ord),cmap,'.',20);
title('Replicates'); 
ylabel(horzcat('Pre-normalization (loss=',num2str(preloss),')'))
legend('location','south');
subplot(2,5,6); 
gscatter(postY(:,1),postY(:,2),MicmaGroupings.data(ord),cmap,'.',20);
ylabel(horzcat('Post-normalization (loss=',num2str(postloss),')'))    
legend('off')
[~,subset,ord]=intersect(Micma.textdata(1,3:end), MicmaClinicalLabels.textdata(2:end,1));
ord=ord+1;
for gid=2:5
    subplot(2,5,gid); 
    gscatter(preY(:,1),preY(:,2),MicmaClinicalLabels.textdata(ord,gid));
    title(MicmaClinicalLabels.textdata{1,gid});
    legend('location','south')
    subplot(2,5,5+gid); 
    gscatter(postY(:,1),postY(:,2),MicmaClinicalLabels.textdata(ord,gid)); 
    legend('off')        
end
suptitle(horzcat('tSNE (perplexity=',num2str(prplx),')'))

