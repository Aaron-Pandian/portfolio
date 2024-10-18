function weightedScore = gradeCalculator(homework, quiz, exam, weight)
    homework = sort(homework,'descend');
    homework = homework(1:end-2,1);
    m = length(homework);
    totalhmkPoints = 0;
    for i = (1:m)
       number = homework(i,1);
       totalhmkPoints = totalhmkPoints + number;
    end
    avgHmk = totalhmkPoints / m;
    
    bonus = quiz(end,1);
    quiz = quiz(1:end-1,1);
    n = length(quiz);
    quiz = sort(quiz,'descend');
    quiz = quiz(1:end-2,1);
    totalPoints = bonus;
    i = 1;
    while i <= length(quiz)
       number = quiz(i,1);
       totalPoints = totalPoints + number;
       i = i + 1;
    end
    avgQuizScore = totalPoints/(n-2)/3*100;
    if avgQuizScore > 100
        avgQuizScore = 100;
    end
    
    midTerm = exam(1,1);
    final = exam(2,1);
    
    hmkWeight = weight(1,1);
    quizWeight = weight(2,1);
    midTermWeight = weight(3,1);
    finalWeight = weight(4,1);
    
    weightedScore = (hmkWeight*avgHmk) + (quizWeight*avgQuizScore) + (midTermWeight*midTerm) + (finalWeight*final);
end
    
    
