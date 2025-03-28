hFig = figure('Color','w');
set(hFig,'unit','normalized','outerPosition',[0 0 1 1])
je = javax.swing.JEditorPane('text/html', '<html><img src="file:/C:\Users\Aritra\OneDrive - Washington University in St. Louis\Lab Workbench\FrankenLab-WashU\Projects\BOVideoProject\BOVideoProjectPrograms\BOVideoTasks\synthetic-02-bees-2x.gif"/></html>');
[hj, hc] =  javacomponent(je,[],hFig);
set(hc, 'unit','normalized','pos',[0.05 0.05 0.2 0.2])