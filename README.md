# My new blog

https://zqfang.github.io/

## How to delopy Hugo cites using github action
see here: [peaceiris](https://github.com/peaceiris/actions-gh-pages)

if repo_name == '<username>.github.io'
do

1. create a source branch and put your source file in this branch
2. use the correct gh-pages.yaml

3. set source
    ```yaml
    name: Publish Site
    on: 
    push:
        branches:
        - source  # default branch
    ```
4. set `publish_branch: master`, or not work
    ```yaml
    - name: Deploy
        uses: peaceiris/actions-gh-pages@v3
        with:
            deploy_key: ${{ secrets.SSH_DEPLOY_KEY }}
            publish_dir: ./public
            publish_branch: master  # deploying branch  
    ```

done