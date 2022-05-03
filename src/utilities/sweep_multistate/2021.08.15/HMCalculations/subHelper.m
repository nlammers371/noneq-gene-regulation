function RNum = subHelper(val,a,subIndices)

    val(subIndices) = val(subIndices)*a;
    valCell = mat2cell(val,size(val,1),ones(1,size(val,2)));
    RNum = RSymFun(valCell{:});