function scc(file1,file2,file3)
	f1=load(file1);
	f2=load(file2);
	G = digraph(f1,f2);
	%p=plot(G);

	%bins = conncomp(G,'Type','strong');

	DG=sparse(f1,f2,true,13748,13748);
	[c, bins] = graphconncomp(DG);

	length(unique(bins));

	fid = fopen(file3,'w+');  % Note the 'wt' for writing in text mode
	for(i=1:length(bins))
		I = find(bins == i);
		fprintf(fid,'%d ',I);  % The format string is applied to each element of a
		fprintf(fid,'\n');
	end
end
