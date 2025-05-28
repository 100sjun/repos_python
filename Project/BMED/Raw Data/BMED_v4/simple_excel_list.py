# 간단한 엑셀 파일 목록 불러오기 예제

import os
import glob

# 방법 1: 가장 간단한 방법
def get_excel_files_simple():
    """현재 폴더의 엑셀 파일 목록을 반환"""
    excel_files = [f for f in os.listdir('.') if f.endswith(('.xlsx', '.xls'))]
    return excel_files

# 방법 2: glob 사용 (추천)
def get_excel_files_glob():
    """glob을 사용하여 엑셀 파일 목록을 반환"""
    return glob.glob('*.xlsx') + glob.glob('*.xls')

# 실행 예제
if __name__ == "__main__":
    # 엑셀 파일 목록 가져오기
    excel_files = get_excel_files_simple()
    
    print(f"현재 폴더에서 {len(excel_files)}개의 엑셀 파일을 찾았습니다:\n")
    
    # 파일 목록 출력
    for i, filename in enumerate(excel_files, 1):
        print(f"{i:2d}. {filename}")
    
    # 리스트로 저장하여 활용
    print(f"\n파일 목록이 'excel_files' 변수에 저장되었습니다.")
    print(f"첫 번째 파일: {excel_files[0] if excel_files else '없음'}")
    print(f"마지막 파일: {excel_files[-1] if excel_files else '없음'}") 